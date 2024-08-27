%Fits kinetic data for rate constants by solving the system of differential
%equations numerically. Handles hidden intermediates and corrects for
%baselines. Input file is data in columns, first x, then y. No headers.
%List species and rate constants separately. Can be named anything, only
%use alphanumeric characters and underscores. Species must be in order of
%data, hidden intermediates must follow their parent. No order needed for
%rate constants. Must set how many species are modeled in each data column.
%Can individually set limits for upper and lower rate constant bounds, and
%initial guesses. Baseline correction is fit value for each set of data
%added to y value to compensate for non-zero baselines. Can individually
%set bounds and initial guesses. Output three files, one with fits and
%data, one with fits with hidden intermediates separated and data, and a
%report file. Numeric solving is much faster than analytical, and testing
%results show the same values are reached.
%First written 2024-01-25
%Updated 2024-02-08: Added hyperlinks to files in output messages
%Updated 2024-03-28: Baseline correction is now subtracted from function
%start point to avoid excess displacement of first point. Non-negative baseline
%lower bounds are now allowed with a warning.
%Updated 2024-05-13: Added capability to fit function starting points.
%Options to enable and set bound distance from initial guess. Also added
%table of Rsquared and RSS values to report file for easier model
%improvement tracking.
%Updated 2024-07-11: Fixed miscount of function start points with hidden 
%intermediates. Implementation of fittable start points failed to account
%for number of species exceeding number of data sets.
clear
close all
tic

[datafilename,datafilepath] = uigetfile('*.txt','Data File');
datafile = fullfile(datafilepath,datafilename);
%datafile = 'Pro7/Rework/Pro7_Kinetic_301_Full.txt';                         %input data. First column x, then y data in columns
basefile = regexprep(datafilename,'.txt','');
outfile = fullfile(datafilepath,basefile);
%outfile = 'Pro7/Rework/Pro7_Kinetic_301_Full_simfit';                 %output base file name. Identifiers and extensions added automatically

%specarray = char('A','B','C');                            %declare each species. One species per equation. Same order as data columns, hidden intermediates follow their parent.
%karray = char('kab','kac');           %declare rate constants. No order needed

%set rate law equations, one for each species. Same order as species.
%{
eqns = {'-kab*A-kac*A';...
        'kab*A';...
        'kac*A';...
        };
%}
%interposit = [1 1 1];                                                   %position array to declare hidden intermediates. Same size as data. List number of equations to fit to each data column
[modelfilename,modelfilepath] = uigetfile('*.txt','Model File');
modelfile = fullfile(modelfilepath,modelfilename);
variables = readcell(modelfile,'Delimiter',' ');

empidx = zeros(size(variables,1),size(variables,2));
emptyloc = zeros(size(variables,1),1);
for i = 1:size(variables,1)
    for j = 1:size(variables,2)
        checkcell = cellfun(@ismissing,variables(i,j),'UniformOutput',false);
        missarray = checkcell{1,end};
        empidx(i,j) = missarray(1,end);
    end
    if sum(empidx(i,:)) == 0
        emptyloc(i,1) = size(variables,2) + 1;
    else
        onepos = find(empidx(i,:)==1);
        emptyloc(i,1) = onepos(1,1);
    end
end

specarray = char(variables{1,1:emptyloc(1,1)-1});
karray = char(variables{2,1:emptyloc(2,1)-1});
eqns = variables(3:size(specarray,1)+2,1);
interposit = variables(size(specarray,1)+3,1);
interposit = str2num(interposit{:}); %#ok<ST2NM>

funcstart = [.09 .33 .17 .09 .04 .05 .04 .02 .03 .05 .01 .03 0];                                                %Initial function points, f(0). One for each equation. If not the same size as #species first data point will be used, hidden intermediates will be zero.
fit_start = 1;                                                                          %1: enable starting point fitting
start_bound = 1;                                                                       %start point bound, fraction above and below initial guess. hard limit is 0 and 1.
lbounds = 0;                                                         %Lower bounds for rate constants. One for each constant. Nonmatching size array will use first value for all.
ubounds = 1000;                                                            %Upper rate constant bound. Same logic as lower
init_guess = .1;                                                   %Initial rate constant guess. One for each rate constant. Missized array uses first value. No entry uses average of bounds
bc_lb = -.01;                                                        %Lower bound for baseline correction. Needs to be negative for best fit. One for each data set. Missized uses first value.                                                        
bc_ub = .1;                                                            %Upper bound. Same logic (doesn't need to be negative)
bc_guess = 0;                                                     %baseline correction first guess. Same logic as bounds. If no values are entered default settings are [-0.01,0.1,0]

%%%%%%
%Beginning of code. Make no changes below this line
%%%%%%

%get sizes of inputs
eqnhold = eqns;                                                             %store raw equations
knum = size(karray,1);                                                      %number of rate constants
specnum = size(specarray,1);                                                %number of species
datnum = size(interposit,2);                                                %number of data sets
%extract data
data = table2array(readtable(datafile));                                    %read data from file
tspan = data(:,1);                                                          %get x axis data
yactual = data(:,2:end);                                                    %get y data

%check number of equations matches number of species
if specnum ~= size(eqns,1)                                                  
    fprintf('Number of equations does not match number of species\n')
    return
end
%check species to fit matches amount of data
if datnum ~= size(yactual,2)
    fprintf('Number of species to fit does not match number of data columns\n')
    return
end
%check position array is correct size
if sum(interposit) ~= specnum
    fprintf('total species set in position array does not match number of species declared\n')
    return
end

%format initial function value
if size(funcstart,2) ~= specnum                                             %check number of initial points equals number of species
    funcstart = zeros(1,specnum);                                           %if not start new array
    funcstart(1,1) = yactual(1,1);                                          %set first value to 1
    for i = 2:datnum
        iter = sum(interposit(1:(i-1)))+1;                                  %get iterator that skips hidden intermediates
        funcstart(iter) = yactual(1,i);                                     %use iterator to get first data point, leaving hidden intermediates at zero
    end
    if fit_start == 1                                                       %check if start point fitting enabled
        f0_lb = funcstart - funcstart.*start_bound;                         %if yes, calculate bounds from initial guess and deviation fraction
        f0_ub = funcstart + funcstart.*start_bound;
        f0_lb(f0_lb < 0) = 0;                                               %put in caps of 0 and 1
        f0_ub(f0_ub > 1) = 1;
        fprintf('fitting function start points\n')
    else
        f0_lb = funcstart;                                                  %no fitting done by setting bounds to guesses
        f0_ub = funcstart;
        fprintf('First data points used as function start points\n')
    end
end
%format rate constant variables
szlb = size(lbounds,2);                                                     %get initial sizes of rate variables
szub = size(ubounds,2);
szig = size(init_guess,2);
if szlb ~= knum                                                             %check lower bounds have correct number of values
    if szlb == 0                                                            %error out if no lower bounds entered
        fprintf('no lower rate bounds entered\n')
        return
    else
        b1 = lbounds(1,1);                                                  %if missized get first value in array
        lbounds = [];                                                       %clear array to avoid oversized array
        lbounds(1,1:knum) = b1;                                             %set correctly sized array to first bound value
    end   
end
if szub ~= size(lbounds,2)                                                  %check if upper bounds match size of correctly sized lower bounds
    if szub == 0
        fprintf('no upper rate bounds entered\n')
        return
    else
        b2 = ubounds(1,1);                                                  %if no, get first bound value
        ubounds = [];                                                       %clear array
        ubounds(1,1:knum) = b2;                                             %resize array and set to first value
    end
end
if szig ~= size(ubounds,2)                                                  %check initial guesses are correct size
    if szig == 0
        bdiff = ubounds-lbounds;                                            %if empty array make guesses average of bounds
        bmid = bdiff./2;
        init_guess = lbounds + bmid;
    else
        i1 = init_guess(1,1);                                               %if missized get first value, clear, and rewrite to correct size
        init_guess = [];
        init_guess(1,1:knum) = i1;
    end
end
boundcheck = lbounds - ubounds;                                             %check upper bounds are greater than lower, get difference
boundcheck(boundcheck < 0) = 0;                                             %larger upper bounds get set to zero
if sum(boundcheck) ~= 0                                                     %if smaller upper bound (non zero) is detected throw error
    fprintf('lower rate bound is greater than upper bound\n')
    return
end
%format baseline correction variables
if size(bc_lb,2) ~= datnum                                                  %check if lower baseline bound is correct size
    if size(bc_lb,2) == 0                                                   %if empty set to -0.01
        bl = -.01;
    else
        bl = bc_lb(1,1);                                                    %if missized set to first value
    end
    bc_lb = [];                                                             %clear array and resize with good values
    bc_lb(1,1:datnum) = bl;
end
if size(bc_ub,2) ~= size(bc_lb,2)                                           %same logic for upper bound
    if size(bc_ub,2) == 0
        bu = .1;                                                            %default value is 0.1
    else
        bu = bc_ub(1,1);
    end
    bc_ub = [];
    bc_ub(1,1:datnum) = bu;
end
if size(bc_guess) ~= size(bc_ub,2)                                          %repeat for initial guess
    if size(bc_guess,2) == 0
        bg = 0;                                                             %default value is 0
    else
        bg = bc_guess(1,1);
    end
    bc_guess = [];
    bc_guess(1,1:datnum) = bg;
end
boundcheck2 = bc_lb - bc_ub;                                                %check upper bound is greater than lower, same as rate bounds
boundcheck2(boundcheck2 < 0) = 0;
if sum(boundcheck2) ~= 0
    fprintf('lower baseline bound is greater than upper bound\n')
    return
end
netcheck = bc_lb;                                                           %check if lower bound is negative. get lower bound
netcheck(netcheck < 0) = 0;                                                 %negative numbers set to zero
if sum(netcheck) ~= 0                                                       %if any nonzero values detected throw error for nonnegative bound
    fprintf('Lower baseline bound contains a nonnegative value.\nIf bad fit occurs try setting lower baseline bound negative\n')
    %return
end

%format equations into function string
for j = 1:size(eqns,1)                                                      %loop through equation set
    workeqn = eqns{j};                                                      %get equation to work on
    for i = 1:size(specarray,1)                                             %loop through species
        specstring = specarray(i,:);
        tarstring = strcat('(?<![\w])',specstring,'(?![\w])');              %find species in equation, using regexp to only get exact matches (species variable surrounded only by non alphanumberic characters)
        repstring = strcat('y(',string(i),')');                             %create substitute string of form y(i). This allows easier variable indexing during fit.
        workeqn = regexprep(workeqn,tarstring,repstring);                   %replace species string with standardized function
    end
    for i = 1:size(karray,1)                                                %loop through rate constants
        kstring = karray(i,:);
        tarstring = strcat('(?<![\w])',kstring,'(?![\w])');                 %same replacement logic as species
        repstring = strcat('k(',string(i),')');                             %standard form is k(i)
        workeqn = regexprep(workeqn,tarstring,repstring);
    end
    eqns{j} = workeqn;                                                      %replace equation with reformatted one.
end
for i=1:size(eqns,1)-1                                                      %add ; to all equations except last. This separates them into individual functions
    eqns{i} = strcat(eqns{i},';');
end
eqnstring = '';
for i=1:size(eqns,1)
    eqnstring = strcat(eqnstring,eqns{i});                                  %concatenate all equations into one string, making system of equations
end
eqnstring = strcat('@(t,y,k,bc)[',eqnstring,']');                           %add function handle, variables, and brackets for proper function formatting               
diffunc = str2func(eqnstring);                                              %convert to proper matlab function

%Fit equations to data
k = optimvar('k',knum,"LowerBound",lbounds,"UpperBound",ubounds);           %set optimization variable for rate constants, bounds define range
bc = optimvar('bc',datnum,"LowerBound",bc_lb,"UpperBound",bc_ub);           %repeat for baseline correction
f0 = optimvar('f0',specnum,'LowerBound',f0_lb,'UpperBound',f0_ub);           %and for function start points
myfcn = fcn2optimexpr(@RtoODE,k,bc,tspan,f0,diffunc,interposit);            %convert function, calling user function to integrate differentials
obj = sum(sum((myfcn - transpose(yactual)).^2));                            %create minimizable function, sum of squares comparison to data
prob = optimproblem("Objective",obj);                                       %create matlab problem to solve
var.k = init_guess;                                                         %format initial guesses into structure
var.bc = bc_guess;
var.f0 = funcstart;
opts=optimoptions(@lsqnonlin,'Display','off','MaxIterations',5000,'OptimalityTolerance',1e-15,'FunctionTolerance',1E-9);        %set solver options
[rsol,sumsq,conditions,output] = solve(prob,var,'Options',opts);            %find best fit rate constants and baseline corrections

%calculate y values with optimized parameters
ysol = transpose(RtoODE(rsol.k,rsol.bc,tspan,rsol.f0,diffunc,interposit));  %evaluate equations with best fit parameters
ysolhid = RtoODEhid(rsol.k,rsol.bc,tspan,rsol.f0,diffunc,interposit);       %evaluate equations, keep hidden intermediates separate.

spcidx = ones(datnum,1);                                                    %preallocate index array
for ii=2:datnum                                                             %find index of species using intermediate position array
    spcidx(ii,1) = sum(interposit(1:ii-1))+1;
end
flist = specarray(spcidx,:);                                                %use index to create list of only plotted species

%calculate rsquared and RSS for each species
rsquare = zeros(1,datnum);
rss = zeros(1,datnum);
rsqreport = strings(datnum,1);                                              %pre-allocate rsqrate arrays
rsqarray = strings(datnum,2);
fitarray = strings(datnum,3);
rssreport = strings(datnum,1);
for i = 1:datnum
    rsquare(i) = Rsquared(yactual(:,i),ysol(:,i));                          %calculate rsquared
    rss(i) = RSS(yactual(:,i),ysol(:,i));                                   %calculate rss
    rsqdisp = sprintf('Rsquared for species %s is %1.4f\n',strtrim(flist(i,:)),rsquare(i));  %format rquared report for species
    rsqarray(i,1) = flist(i,:);                                             %put species name in array
    rsqarray(i,2) = rsquare(i);                                             %put rsquared value in array
    rsqreport(i) = string(rsqdisp);                                         %store rsquared report
    rssreport(i) = string(sprintf('RSS for species %s is %1.6f\n',strtrim(flist(i,:)),rss(i)));         %format rss report for each species
    fitarray(i,1) = flist(i,:);                                             %form fitting report array   
    fitarray(i,2) = rsquare(i);
    fitarray(i,3) = rss(i);
end
rln1 = 'R squared';                                                         %title
rln2 = sprintf('\n%6s %10s',transpose(rsqarray));                           %form array of rquared values for each species
rsqlist = append(rln1,rln2);                                                %attach title to array
disp(rsqlist)                                                               %display rsquared values
infostore = cell(datnum,1);                                                 %preallocate report array

%calculate overall Rsquared
Re = yactual(:,1:end)-ysol;                                                 %Rsquared calculation for comparison to previous data
Resq = Re.^2;
Sres = sum(Resq,1);
Ryhat = sum(yactual(:,1:end),1)./size(yactual,1);
Ryf = (yactual(:,1:end)-Ryhat).^2;
Stot = sum(Ryf,1);
StotT = sum(Stot);                                                          %sum residual average for overall fit
SresT = sum(Sres);                                                          %sum RSS for all fits to judge overall fit
RsquareT = 1-(SresT/StotT);                                                 %overall Rsquared
rsqtot = sprintf('Overall Rsquared is %1.6f\n',RsquareT);                   %format and store values for overall Rsq and RSS
rsstot = sprintf('Overall RSS is %3.6f\n',SresT);
disp(rsqtot)                                                                %display overall Rsq
infostore{end+1} = rsqtot;                                                  %write rsq and rss to report array
infostore{end+1} = rsstot;

%Form table of Rsquare and RSS
fitarray(end+1,:) = 'Overall';                                              %add last row for overall fit
fitarray(end,2) = RsquareT;
fitarray(end,3) = SresT;
fln1 = 'R squared and RSS';                                                 %header
fln2 = sprintf('\n%10s %15s %15s',transpose(fitarray));                     %format data array
fitreport = append(fln1,fln2);                                              %combine header and data array
infostore{end+1} = fitreport;                                               %store in report cell

%Print value of each k value
klistarray = strings(knum+1,2);                                             %pre-allocate rate constant array
for ii=1:knum
    klistarray(ii+1,1) = karray(ii,:);                                      %pull name of each rate constant
    klistarray(ii+1,2) = rsol.k(ii);                                        %pull value of each rate constant
end
kln1 = 'Rate Constants';                                                    %title
kln2 = sprintf('%6s %20s\n',transpose(klistarray));                         %format rate constant value with name
kreport = append(kln1,kln2);                                                %form report array
disp(kreport)



%Form report array
for ii=1:datnum                                                             %loop through species
    rsqout = rsqreport(ii);                                                 %for each species add rsq...
    rssout = rssreport(ii);                                                 %rss...
    baseout = sprintf('baseline correction for species %s is %f\n',strtrim(flist(ii,:)),rsol.bc(ii));   %baseline correction...
    %eqnout = sprintf('%s\n',eqns{ii});                                     %and equation
    report = append(rsqout,rssout,baseout);                                 %combine together into one array
    infostore{ii} = report;                                                 %store for each species
end

%form fit report
initstore = strings(specnum,1);
for i = 1:specnum                                                           %loop through species, and format string for function start point
    initstore(i,:) = sprintf('function %s starting point: %f\n',strtrim(specarray(i,:)),rsol.f0(i,1));
end
infostore{end+1} = kreport;                                                 %add rate constants...
infostore{end+1} = sprintf('Input data file: %s\n',datafile);               %input data file...
infostore{end+1} = initstore;                                               %starting points...
infostore{end+1} = specarray;                                               %species...
infostore{end+1} = karray;                                                  %rate constants
eqnstore = strings(specnum,1);                                              
for i = 1:specnum
    eqnstore(i,:) = eqnhold{i};                                             %write each input equation to array
end
infostore{end+1} = eqnstore;                                                %and report input equations

%Plot data and fits
colorlist = distinguishable_colors(datnum);                                 %define color list to use in plots
figure(1)
g1 = zeros(1,datnum);
for ii = 1:datnum                                                           %plot y data by looping
    g1(ii) = plot(tspan,yactual(:,ii),'.','MarkerSize',20,'Color',colorlist(ii,:));
    hold on
end
for ii = 1:datnum                                                           %plot fits, using same color list as data to match
    g2 = plot(tspan,ysol(:,ii),'LineStyle','-','Color',colorlist(ii,:),'LineWidth',2);
end
titlestring = sprintf('Overall Data Fit (R^{2} = %1.4f)',RsquareT);
title(titlestring)                                                   %set title
legend(g1,flist);                                                           %set plot legend only for data
xlabel('Time');                                                             %set plot axis labels
ylabel('Abundance');
hold off
%{
%Plot hidden intermediates with data
excolorlist = zeros(specnum,1);                                             %preallocate color array
for ii = 1:datnum                                                           %form list to match color of each component equation to corresponding data
    specount = sum(interposit(1:ii));                                       %form index value by looping through intermediate declaration array and taking the sum at each index
    excolorlist(specount,1) = ii;                                           %count up to number of species, skipping an index for each component of a species
end
for ii = 2:specnum                                                          %loop back through number of component equations to backfill skipped indices
    ri = specnum+1-ii;                                                      %create iterator to work through array in reverse
    if excolorlist(ri,1) == 0                                               %check if value was skipped
        excolorlist(ri,1) = excolorlist(ri+1,1);                            %if yes, fill with value that follows it in the array. This assigns the same value for all intermediates within a species
    end                                                                     %this results in array where each component of a species has the same index number, and so can call the same color
end
figure(2)
for ii = 1:specnum                                                          %plot fit data by looping through, picking colors with color index created above
    g3 = plot(tspan,ysolhid(:,ii),'LineStyle','--','Color',colorlist(excolorlist(ii),:),'LineWidth',2);
    hold on
end
g4 = zeros(1,datnum);
for ii = 1:datnum                                                           %plot y data by looping through it, picking colors from list.
    g4(ii) = plot(tspan,yactual(:,ii),'.','MarkerSize',20,'Color',colorlist(ii,:));
end
legend(g4,flist);                                                           %set legend only for data
title('Fitted Data including Hidden Intermediates')                         %set title, x and y labels for plot
xlabel('Time');
ylabel('Abundance');
hold off
%}
%Write data to file
outdata = zeros(size(yactual,1),2*datnum+1);                                %set up output data file
outdata(1:end,1) = tspan;                                                   %add in x data...
outdata(1:end,2:datnum+1) = yactual;                                        %y data...
outdata(1:end,(datnum+2):end) = ysol;                                       %and fit data
dataout = strcat(outfile,'_fit.txt');                                       %form output file name from base output name
writematrix(outdata,dataout)                                                %write data to file
fprintf('Wrote output data to %s\n',dataout);                               %display file location

%Write hidden intermediate data to file
outdata2 = zeros(size(tspan,1),datnum+specnum+1);                           %set up output data file
outdata2(:,1) = tspan;                                                      %add in x data...
outdata2(:,2:datnum+1) = yactual;                                           %y data...
outdata2(:,(datnum+2):end) = ysolhid;                                       %and fits of individual equations
dataout2 = strcat(outfile,'_indiv_data.txt');                               %form output file name from base output name
%writematrix(outdata2,dataout2)                                              %write data to file
%fprintf('Wrote individual fit data to <a href="matlab: winopen(dataout2)">%s</a>\n',dataout2);                      %display file location

%Write report to file
reportout = strcat(outfile,'_report.txt');                                  %form report file name from base output name
%writecell(infostore,reportout,'Delimiter',' ');                             %write report to file
%fprintf('Wrote fitting report to <a href="matlab: winopen(reportout)">%s</a>\n',reportout);                          %display file location

toc
function R = Rsquared(y,yf)                                                 %Calculates R-squared from signal (y) and fit (yf)
    e = y-yf;                   
    esq = e.^2;
    esq(isnan(esq))=0;                                                      %NaN values occur from blank spreadsheet inputs or matlab rounding tiny numbers to 0
    y(isnan(y))=0;                                                          %Setting them to zero here allows calculation to proceed, otherwise it will not converge
    yhat = sum(y)/size(y,1);
    yf = (y-yhat).^2;
    Stot = sum(yf);
    Sres = sum(esq);
    R = 1-(Sres/Stot);
end
function hidata = hidintdata(y,positionarray)                               %Adds hidden intermediates to their parent data
    datacol = size(positionarray,2);
    hidata = zeros(size(y,1),datacol);
    for i = 1:datacol
        if i==1                                                             %set indices of equations for first species
            ledx = 1;
            hedx = positionarray(1);
        else                                                                %use intermediate position array to index hidden intermediates with their species
            ledx = sum(positionarray(1:(i-1)))+1;
            hedx = sum(positionarray(1:i));
        end
        hidata(:,i) = sum(y(:,ledx:hedx),2);
    end
end
function solpts = RtoODE(k,bc,tspan,y0,diffunc,positionarray)               %Numerically integrates system of differential equations, calcuates y values, and adds a baseline
    idx = hididx(positionarray);                                            %get index of real species (excludes hidden intermediates)
    y0(idx,1) = y0(idx,1) - bc;                                             %subtract baseline from starting point (y0 usually first data point which inherently includes baseline)
    sol = ode45(@(t,y)diffunc(t,y,k),tspan,y0);                             %numerically solve system of differential equations
    solptsraw = deval(sol,tspan);                                           %evaluate solution at each x point
    solpts1 = hidintdata(transpose(solptsraw),positionarray);               %combine hidden intermediate fits
    solpts = transpose(solpts1) + bc;                                       %add baseline correction
end
function solpts = RtoODEhid(k,bc,tspan,y0,diffunc,positionarray)            %Numerically integrates system of differentials, leaving hidden data expressed
    idx = hididx(positionarray);                                            %get index of real species (excludes hidden intermediates)
    y0(idx,1) = y0(idx,1) - bc;                                  %subtract baseline from starting point (y0 usually first data point which inherently includes baseline)
    sol = ode45(@(t,y)diffunc(t,y,k),tspan,y0);                             %solve system of differential equations
    solptsraw = deval(sol,tspan);                                           %evaluate solution at each x point
    bcindiv = zeros(1,size(solptsraw,1));                                   %set baseline array for all fits
    bcindiv(1,1) = bc(1);                                                   %first slot is first baseline value
    for i = 2:size(positionarray,2)
        iter = sum(positionarray(1:(i-1)))+1;                               %loop through baseline array, leaving hidden intermediates at zero, setting baseline for parent
        bcindiv(iter) = bc(i);
    end
    solpts = transpose(solptsraw) + bcindiv;                                %add baseline to fits
end
function ledx = hididx(positionarray)                                       %creates index of real species in total species list
    datacol = size(positionarray,2);                                        %get number of species
    ledx = ones(1,datacol);                                                 %pre-allocate index with ones, gives starting index
    for i = 2:datacol
        ledx(:,i) = sum(positionarray(1:(i-1)))+1;                          %use intermediate position array to create index list of real species
    end
end
function Q = RSS(y,yf)                                                      %Calculates RSS    
    e = y-yf;                   
    esq = e.^2;
    esq(isnan(esq))=0;                  %NaN values occur from blank spreadsheet inputs or matlab rounding tiny numbers to 0
    %y(isnan(y))=0;                     %Setting them to zero here allows calculation to proceed, otherwise it will not converge
    Sres = sum(esq);
    Q = Sres;
end
function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);
% Copyright 2010-2011 by Timothy E. Holy
% https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end
  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
end
function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
end
function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1
    c = rgbspec(k,:);
  elseif length(c)>2
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
end
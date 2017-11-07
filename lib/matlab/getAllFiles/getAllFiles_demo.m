%% _getAllFiles_ demo
% The _getAllFiles_ function recursively collects a list of files from a
% folder tree, allowing you to specify the selection criteria for which
% files are collected and how they are formatted. For these examples, we'll
% be using the main MATLAB toolbox path.
%%

rootPath = 'C:\Program Files\MATLAB\R2016b\toolbox\matlab';
format compact

%% The 'FileFilter' option
% We can specify a
% <https://www.mathworks.com/help/matlab/matlab_prog/regular-expressions.html
% regular-expression> pattern to filter the file names on, collecting a
% list of files that match. Here's an example that uses the 'FileFilter'
% option to recursively find every file with a '.m' extension:

fileList = getAllFiles(rootPath, 'FileFilter', '\.m$');
fprintf('%d files found.\n', size(fileList, 1));
fprintf('%s\n', fileList{1:5}, '...');

%%
% It's a pretty long list, so I've only shown the first 5 files it finds.
% Notice they are listed with the full path prepended by default.

%%
% If you have multiple match expressions to filter on, you can use
% <https://www.mathworks.com/help/matlab/matlab_prog/regular-expressions.html#f0-42884
% grouping operators> to include them all in one expression. For example,
% this will find every '.jpg', '.bmp', and '.tif' file:

fileList = getAllFiles(rootPath, 'FileFilter', '\.(jpg|bmp|tif)$');
fprintf('%d files found.\n', size(fileList, 1));
fprintf('%s\n', fileList{1:5}, '...');

%% The 'Depth' option
% If we don't want to search quite so far down the folder tree, we can
% limit the search depth with the 'Depth' option. Let's see how many '.m'
% files are in the root folder:

fileList = getAllFiles(rootPath, 'FileFilter', '\.m$', 'Depth', 0);
fprintf('%d files found.\n', size(fileList, 1));

%%
% Looks like none are. They are all contained in subfolders. Let's see how
% many are located in just the immediate subfolders of the root folder:

fileList = getAllFiles(rootPath, 'FileFilter', '\.m$', 'Depth', 1);
fprintf('%d files found.\n', size(fileList, 1));

%% The 'PrependPath' option
% Maybe we just want the file names, but don't care about the absolute
% paths. In this case, we just set the 'PrependPath' option to _false_:

fileList = getAllFiles(rootPath, 'FileFilter', '\.m$', ...
                                 'PrependPath', false);
fprintf('%d files found.\n', size(fileList, 1));
fprintf('%s\n', fileList{1:5}, '...');

%% The 'ValidateFcn' option
% Sometimes we might want to select files based on a more complicated
% criteria than just what's in their names. In this case, we can use the
% 'ValidateFcn' option to specify a function that is to be run on each file
% found. This function should accept as input a structure of the form
% returned by the |dir| function and return a logical value (|true| to
% collect it in the list, |false| to ignore it). First, let's find all the
% '.png' files:

fileList = getAllFiles(rootPath, 'FileFilter', '\.png$');
fprintf('%d files found.\n', size(fileList, 1));

%%
% Now, we can specify an anonymous function that gets the byte size of each
% file and returns |true| for only those greater than 250KB:

bigFcn = @(s) (s.bytes > 512^2);
fileList = getAllFiles(rootPath, 'FileFilter', '\.png$', ...
                                 'ValidateFcn', bigFcn);
fprintf('%s\n', fileList{:});

%%
% Just the one.
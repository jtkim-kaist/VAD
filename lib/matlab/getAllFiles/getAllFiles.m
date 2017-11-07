function fileList = getAllFiles(rootPath, varargin)
%getAllFiles   Recursively collect files within a folder.
%   FILELIST = getAllFiles(ROOTPATH) will search recursively through the
%   folder tree beneath ROOTPATH and collect a cell array FILELIST of all
%   files it finds. The list will contain the absolute paths to each file
%   starting at ROOTPATH.
%
%   FILELIST = getAllFiles(ROOTPATH, 'PropertyName', PropertyValue, ...)
%   will modify how files are selected based on the property/value pairs
%   specified. Valid properties that the user can set are:
%
%     'Depth'       - A non-negative integer value for the maximum folder
%                     tree depth that getAllFiles will search through. A
%                     value of 0 will only search in ROOTPATH, a value of 1
%                     will search in ROOTPATH and its subfolders, etc..
%                     Default (and maximum allowable) value is the current
%                     recursion limit set on the root object (i.e.
%                     get(0, 'RecursionLimit')).
%     'FileFilter'  - A string defining a regular-expression pattern that
%                     will be applied to the file name. Only files matching
%                     the pattern will be returned in FILELIST. Default is
%                     '' (i.e. all files are returned).
%     'PrependPath' - A logical value determining if the full path from
%                     ROOTPATH to the file is prepended to each file in
%                     FILELIST. The default TRUE will prepend the full
%                     path, otherwise just the file name is returned.
%     'ValidateFcn' - A handle to a function that takes as input a
%                     structure of the form returned by the DIR function
%                     and returns a logical value. This function will be
%                     applied to all files found and only files that have
%                     a TRUE return value will be returned in FILELIST.
%                     Default is [] (i.e. no validation done).
%
%   Examples:
%
%     1) Find all '.m' files:
%
%        fileList = getAllFiles(rootPath, 'FileFilter', '\.m$');
%
%     2) Find all '.jpg', '.png', and '.tif' files:
%
%        fileList = getAllFiles(rootPath, ...
%                               'FileFilter', '\.(jpg|png|tif)$');
%
%     3) Find all '.m' files in the given folder and its subfolders:
%
%        fileList = getAllFiles(rootPath, 'Depth', 1, ...
%                                         'FileFilter', '\.m$');
%
%     4) Find all '.m' files, returning only the file names:
%
%        fileList = getAllFiles(rootPath, 'FileFilter', '\.m$', ...
%                                         'PrependPath', false);
%
%     5) Find all '.jpg' files with a size of more than 1MB:
%
%        bigFcn = @(s) (s.bytes > 1024^2);
%        fileList = getAllFiles(rootPath, 'FileFilter', '\.jpg$', ...
%                                         'ValidateFcn', bigFcn);
%
%   See also dir, regexp.

% Author: Ken Eaton
% Version: MATLAB R2016b
% Last modified: 3/6/17
% Copyright 2017 by Kenneth P. Eaton
%--------------------------------------------------------------------------

  % Create input parser (only have to do this once, hence the use of a
  %   persistent variable):

  persistent parser
  if isempty(parser)
    recursionLimit = get(0, 'RecursionLimit');
    parser = inputParser();
    parser.FunctionName = 'getAllFiles';
    parser.PartialMatching = true;
    addRequired(parser, 'rootPath', ...
                @(s) validateattributes(s, {'char'}, {'nonempty'}));
    addParameter(parser, 'Depth', recursionLimit, ...
                 @(s) validateattributes(s, {'numeric'}, ...
                                         {'scalar', 'nonnegative', ...
                                          'nonnan', 'integer', ...
                                          '<=', recursionLimit}));
    addParameter(parser, 'FileFilter', '', ...
                 @(s) validateattributes(s, {'char'}, {'row'}));
    addParameter(parser, 'PrependPath', true, ...
                 @(b) validateattributes(b, {'logical'}, {'scalar'}));
    addParameter(parser, 'ValidateFcn', [], ...
                 @(f) validateattributes(f, {'function_handle'}, ...
                                         {'scalar'}));
  end

  % Parse input and recursively find files:

  parse(parser, rootPath, varargin{:});
  fileList = getAllFiles_core(parser.Results.rootPath, 0, ...
                              rmfield(parser.Results, 'rootPath'));

end

%~~~Begin local functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%--------------------------------------------------------------------------
function fileList = getAllFiles_core(rootPath, depth, optionStruct)
%
%   Core recursive function to find files.
%
%--------------------------------------------------------------------------

  % Find all valid subdirectories in current directory:

  dirData = dir(rootPath);
  dirIndex = [dirData.isdir];
  subDirs = {dirData(dirIndex).name};
  subDirs = subDirs(~ismember(subDirs, {'.', '..'}));
  subDirs = fullfile(rootPath, subDirs);
  nSubs = numel(subDirs);

  % Find all files in the current directory:

  dirData = dirData(~dirIndex);
  fileList = {dirData.name}.';

  % Perform culling and processing on file list:

  if ~isempty(fileList)

    % Apply file name filter, if specified:

    fileFilter = optionStruct.FileFilter;
    if ~isempty(fileFilter)
      filterIndex = ~cellfun(@isempty, regexp(fileList, fileFilter));
      dirData = dirData(filterIndex);
      fileList = fileList(filterIndex);
    end

    if ~isempty(fileList)

      % Apply validation function to file list, if specified:

      validateFcn = optionStruct.ValidateFcn;
      if ~isempty(validateFcn)
        fileList = fileList(arrayfun(validateFcn, dirData));
      end

      % Prepend full path to file list, if specified:

      if (~isempty(fileList) && optionStruct.PrependPath)
        fileList = fullfile(rootPath, fileList);
      end

    end

  end

  % Collect file lists from subdirectories:

  if (depth < optionStruct.Depth)
    fileList = {fileList; cell(nSubs, 1)};
    for iSub = 1:nSubs
      fileList{iSub+1} = getAllFiles_core(subDirs{iSub}, depth+1, ...
                                          optionStruct);
    end
    fileList = vertcat(fileList{:});
  end

end

%~~~End local functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
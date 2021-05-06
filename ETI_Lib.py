import os
import sys
import time
import datetime
import logging
import subprocess
# To parse and use setting file
from configparser import ConfigParser

# Biopython
from Bio import Entrez
from Bio.Entrez import Parser as ps


def getBaseName(path):
    """The program filename without extension"""
    return os.path.basename(path).split('.')[0]


def beingOrExit(prgPath):
    """Check if the given path exists and eventually halt the execution"""
    if not os.path.exists(prgPath):
        print('The file ', prgPath, 'does not exist, check the path')
        sys.exit()


def run(prgName):
    process = subprocess.run([prgName],
                             shell=True,
                             check=False,
                             capture_output=True,
                             text=True)
    return process


def pathFilename(path, filename):
    """Join path and filename"""
    return os.path.join(path, filename)


def setPath(path):
    """Check the existence of a path and build it if necessary."""
    print('setPath path', path)
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def pathJoinCheck(dir2add, rootPath='.'):
    """Join the dir2add to the rootpath. It builds the path if necessary."""
    path = os.path.join(rootPath, dir2add)
    setPath(path)
    return path


def setWorkingDirs(wrkDir=None, dizDirNames=None):
    """Set the working directories as recived from input dictionary.

    If arguments are omitted, current working directory is set as
    workingDirectory and no directory is nested.
    If wrkDir does not exists, the script create it.
    wrkDir/dir1
          /dir2
          /dir3

    Keyword arguments:
     wrkDir: The working directory path. if None it is set to cwd
     dizDirNames: The dictionary containing the names of the nested
                  directories, in the form: {'dirName', 'dirPath'}

    Return:
    (wrkDirPath, dizDirPaths) where:
    wrkDirPath: The working dir path
    dizDirPaths: The dictionary containing the paths of the nested directories,
     in the form: {'dirName', 'dirPath'}
    """
    cwDir = os.getcwd()
    dizDirPaths = {}
    if wrkDir is None:
        wrkDirPath = cwDir + os.sep
    else:
        if not os.path.exists(wrkDir):
            print('the working directory', wrkDir, 'does not exist. Fix it.')
            sys.exit()
        else:
            wrkDirPath = os.path.abspath(wrkDir)
            if wrkDirPath[-1] != os.sep:
                wrkDirPath += os.sep
    if dizDirNames is not None:
        for el in dizDirNames:
            dizDirPaths[el] = pathJoinCheck(dizDirNames[el], wrkDirPath)
    return (wrkDirPath, dizDirPaths)


def dirPathFromConfig(configObj=None,
                      configParamName=None,
                      default=None,
                      wrkDir=None):
    if wrkDir is None:
        wrkDir = os.getcwd()
    dir = configObj.get('parameters', configParamName, fallback=None)
    print('dir', dir)
    if dir is None:
        dirPath = pathJoinCheck(default, wrkDir)
    else:
        print('absPath', os.path.abspath(dir))
        if dir[0] == '~':
            dir = os.path.join(os.path.expanduser(dir[0]), dir[2:])
        print('dir', dir)
        print('absPath', os.path.abspath(dir))

        if not os.path.exists(dir):
            print('the working directory', dir, 'does not exist. Fixing it.')
            dirPath = setPath(dir)
        else:
            print('the path exists')
            dirPath = os.path.abspath(dir)
            if dirPath[-1] != os.sep:
                dirPath += os.sep
    return dirPath


def initFileSystem(configObj):
    """Set the working directories as recived from input dictionary.

    If arguments are omitted, current working directory is set as
    workingDirectory and no directory is nested.
    If wrkDir does not exists, the script create it.
    wrkDir/dir1
          /dir2
          /dir3

    Keyword arguments:
     wrkDir: The working directory path. if None it is set to cwd
     dizDirNames: The dictionary containing the names of the nested
                  directories, in the form: {'dirName', 'dirPath'}

    Return:
    (wrkDirPath, dizDirPaths) where:
    wrkDirPath: The working dir path
    dizDirPaths: The dictionary containing the paths of the nested directories,
     in the form: {'dirName', 'dirPath'}
    """
    cwDir = os.getcwd()
    dizDirPaths = {}
    dizDirPaths['wrk'] = dirPathFromConfig(configObj,
                                           'working_dir_path',
                                           default=cwDir,
                                           wrkDir=None)
    dizDirPaths['ete'] = dirPathFromConfig(configObj,
                                           'ete_path',
                                           default='Ete',
                                           wrkDir=dizDirPaths['wrk'])
    dizDirPaths['down'] = dirPathFromConfig(configObj,
                                            'download_path',
                                            default='Downloads',
                                            wrkDir=dizDirPaths['wrk'])
    dizDirPaths['log'] = dirPathFromConfig(configObj,
                                           'log_folder',
                                           default='Logs',
                                           wrkDir=dizDirPaths['wrk'])
    return dizDirPaths


def getTimeStamp():
    return time.strftime('%Y%m%d%H%M%S', time.localtime())


def logFileOpen(logDIR=None, timeStamp=None, aim=None, basename=None):
    if logDIR is None:
        logDIR = os.getcwd() + os.sep
    if timeStamp is None:
        timeStamp = getTimeStamp()
    if aim is None:
        aimStr = ''
    else:
        aimStr = '_' + aim
    if basename is None:
        basenameStr = ''
    else:
        basenameStr = basename + '_'
    logFileName = os.path.join(logDIR,
                               basenameStr + timeStamp + aimStr + '.log')
    logStream = open(logFileName, mode='w')
    return logStream


def toLog(logStream, string):
    logStream.write(string + '\n')


def logInit(logFolder=None, timeStamp=None):
    """Initialize the logging activity"""
    #  try:
    #      # In case the log folder in the tool folder doesn't exist
    #      if not os.path.exists("./Log/"):
    #          # It will try to create it
    #          os.makedirs("./Log/")
    #          # In case everything goes well, log file will be created in the log folder
    #          # The name of the file will be the date and the time at which the tool has started
    #      log_folder = './{0}/{1}.{2}'.format(
    #          "Log", str(datetime.datetime.now().strftime("%d-%m-%y_at_%H-%M")),
    #          "log")
    #
    #  except OSError:
    #      # In case of denied permissions or other cases where we can't create the log folder, the file will be in the same
    #      # folder of the tool
    #      # The name of the file will be the date and the time at which the tool has started
    #      log_folder = '{0}.{1}'.format(
    #          str(datetime.datetime.now().strftime("%d-%m-%y_at_%H-%M")), "log")

    # Creating the log file
    logging.basicConfig(filename=pathFilename(logFolder,
                                              'extaxsi_' + timeStamp + '.log'),
                        format='%(asctime)s  %(levelname)s -> %(message)s',
                        level=logging.DEBUG)


# DA CAMBIARE PRIMA DELLA PUBBLICAZIONE (DEBUG -> WARNING) !!!!!!!!!!!!!!!!!


# Creates the path where the files will be downloaded
def create_folder(folder_path):
    logging.info(' Searching or creating the download path')
    try:
        if not os.path.exists(
                folder_path):  # If the folder doesn't exist it will make it
            os.makedirs(folder_path)
            logging.info(' Path created!')
            return True
        else:
            logging.info(' Path already exist')
            return True
    except OSError:
        print('Error: Creating directory ' + folder_path,
              '\nFiles will be saved in the tool directory')
        logging.warning(
            ' Error at creating directory << %s >> try to change your path in the settings file'
            % folder_path)
        return False


# Rename files name
def rename_file(folder_path, rename_me, extension):

    while "/" in rename_me:

        if "/" in rename_me:
            rename_me = rename_me[rename_me.find("/") + 1:]

    if not os.path.exists(folder_path):
        rename_me = '{0}{1}'.format(
            rename_me.replace(" ", "_").replace("(", "").replace(")", ""),
            extension)

    else:
        rename_me = '{0}{1}{2}'.format(
            folder_path,
            rename_me.replace(" ", "_").replace("(", "").replace(")", ""),
            extension)

    return rename_me


# Initializing the Filesystem:
runningTimeStamp = getTimeStamp()

# Reading the config file
config = ConfigParser()
config.read('settings.ini')

dWrkDirs = initFileSystem(config)

logInit(logFolder=dWrkDirs['log'], timeStamp=runningTimeStamp)

# Here goes your email so NCBI can contact you in case of necessity
logging.info(" Getting entrez_mail parameter from config file")
Entrez.email = config.get('parameters', 'entrez_email')

if Entrez.email == "example@example.com":
    print("Please enter your email in the settings.ini file")
    logging.info(" entrez_mail parameter is not defined (default value)")

else:
    logging.info(" Getting entrez_mail parameter completed")

# Here goes your api key so
logging.info(" Getting api_key parameter from config file")
if config.get('parameters', 'api_key') == "none":
    logging.info(
        " No api_key parameter entered, for best results enter your ncbi api key in the settings.ini file"
    )

else:
    Entrez.api_key = config.get('parameters', 'api_key')
    logging.info(" Getting api_key parameter completed")

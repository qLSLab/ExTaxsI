import os
import datetime
import logging


def logInit():
    """Initialize the logging activity"""
    try:
        # In case the log folder in the tool folder doesn't exist
        if not os.path.exists("./Log/"):
            # It will try to create it
            os.makedirs("./Log/")
            # In case everything goes well, log file will be created in the log folder
            # The name of the file will be the date and the time at which the tool has started
        log_folder = './{0}/{1}.{2}'.format(
            "Log", str(datetime.datetime.now().strftime("%d-%m-%y_at_%H-%M")),
            "log")

    except OSError:
        # In case of denied permissions or other cases where we can't create the log folder, the file will be in the same
        # folder of the tool
        # The name of the file will be the date and the time at which the tool has started
        log_folder = '{0}.{1}'.format(
            str(datetime.datetime.now().strftime("%d-%m-%y_at_%H-%M")), "log")

    # Creating the log file
    logging.basicConfig(filename=log_folder,
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

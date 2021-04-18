import sys
import os
import time

# To check which OS we are using
import platform

# To give tab functionality while searching for files
# in bash (Not working in windows)
if platform.system() != "Windows":
    import readline

import ETI_Lib_DB as ETIDB
import ETI_Lib_Statistics as ETILSTAT

class color:
    # Colors
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    PURPLE = '\033[95m'
    CYAN = '\033[96m'

    # Text format
    BOLD = '\033[1m'
    ITALIC = '\33[3m'
    UNDERLINE = '\033[4m'
    BLINK = '\33[5m'

    # Background colors
    BLACKBG = '\33[40m'
    REDBG = '\33[41m'
    GREENBG = '\33[42m'
    YELLOWBG = '\33[43m'
    BLUEBG = '\33[44m'
    VIOLETBG = '\33[45m'
    BEIGEBG = '\33[46m'
    WHITEBG = '\33[47m'

    # End of format\colors
    END = '\033[0m'


def DNA_drawing():
    print("\n")
    print(color.GREEN + "      .%%%                " + color.RED + "*%%(  " +
          color.END)
    print(color.GREEN + "       %%%#  %%%%%%" + color.RED + "%%%%%% .%%%  " +
          color.END)
    print(color.GREEN + "        %%%%              " + color.RED + " %%%/ " +
          color.END)
    print(color.GREEN + "         #%%%#  %%%%%" + color.RED + "%%%%% %%%/ " +
          color.END)
    print(color.GREEN + "           %%%%%           " + color.RED + "%%%. " +
          color.END)
    print(color.GREEN + "             %%%%%.  %%" + color.RED + "%% %%%%  " +
          color.END)
    print(color.GREEN + "               #%%%%/   " + color.RED + " %%%%   " +
          color.END)
    print(color.GREEN + "                 *%%%# " + color.RED + ",%%%%    " +
          color.END)
    print(color.GREEN + "                     %" + color.RED + "%%%%*     " +
          color.END)
    print(color.RED + "                   *%%%% " + color.GREEN +
          "%*         " + color.CYAN +
          "-------------------------------------------------------" +
          color.END)
    print(color.RED + "                 %%%%%* " + color.GREEN +
          "*%%%        " + color.CYAN + "----------------- WELCOME TO " +
          color.BOLD + color.UNDERLINE + "ExTaxsI" + color.END + color.CYAN +
          "! -----------------" + color.END)
    print(color.RED + "              .%%%%%     " + color.GREEN +
          "#%%%       " + color.CYAN +
          "-------------------------------------------------------" +
          color.END)
    print(color.RED + "            *%%%%#  *%%" + color.GREEN +
          "%  #%%%      " + color.CYAN +
          "--------- Exploration Taxonomies Informations ---------" +
          color.END)
    print(color.RED + "          ,%%%%*           " + color.GREEN +
          "%%%,     " + color.CYAN +
          "---- Tool to build and visualize amplicon database ----" +
          color.END)
    print(color.RED + "         %%%%*  %%%%%" + color.GREEN +
          "%%%%% %%%(     " + color.CYAN +
          "-------------------------------------------------------" +
          color.END)
    print(color.RED + "       .%%%#              " + color.GREEN + ".&%%*   " +
          color.END)
    print(color.RED + "      .%%%/  %%%%%" + color.GREEN + "%%%%%% ./%%%    " +
          color.END)
    print(color.RED + "      %%%#              " + color.GREEN + "..*%%*    " +
          color.END)
    print(color.RED + "      %%%  %%%%%" + color.GREEN + "%%%%% ..#%%%(     " +
          color.END)
    print(color.RED + "     .%%%             " + color.GREEN + "*%%%%       " +
          color.END)
    print(color.RED + "      %%%,  %%%" + color.GREEN + "%%%  (%%%%,        " +
          color.END)
    print(color.RED + "      (%%%        " + color.GREEN + "%%%%%.          " +
          color.END)
    print(color.RED + "       #%%% %" + color.GREEN + "% ,%%%%%             " +
          color.END)
    print(color.RED + "        ,%%, " + color.GREEN + " (%%%(               " +
          color.END)
    print(color.RED + "           *" + color.GREEN + "/%%%,                 " +
          color.END)
    print(color.GREEN + "          .%%%%.                              " +
          color.END)



def update_progress(work_needed, work_done):
    """Function to animate a progress bar

    Keyword arguments:
        work_needed
        wor_done


    Return:
         None
    """

    # scritto da adam non funziona:
    # if not is_integer(work_done) or work_done <= 0:
    #    return

    bar_len = 40  # Modify this to change the length of the progress bar
    status = ""  # Initializing the string that will give us the status of the progression
    progress = work_needed / work_done

    if isinstance(progress, int):
        progress = float(progress)

    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"

    if progress < 0:  # progress can't be negative
        progress = 0
        status = "Halt...\r\n"

    if progress >= 1:
        progress = 1
        status = "Done...\r\n"

    block = int(round(bar_len * progress))  # % value of the progress
    text = '\r{6}Loading:{7} {5}[{0}]{7} {1}% || {2}/{3} {6}{4}{7}'.format(
        '█' * block + '~' * (bar_len - block), round(progress * 100, 1),
        work_needed, work_done, status, color.GREEN, color.BOLD, color.END)
    sys.stdout.write(text)
    sys.stdout.flush()


def clear():
    if platform.system() == "Windows":
        os.system('cls')
    else:
        os.system('tput reset')


def main_menu():
    DNA_drawing()
    print("\n Which module do you want to use? \n")
    print(" 1 - Database creation module: taxonomy and fasta files download ")
    print(
        " 2 - Statistical module: Scatter plot and world map from taxonomy files or queries "
    )
    # print(" 3 - BLAST functions: create or use a profile and test your new database ")
    print(
        " 3 - Taxonomy ID converter: convert a single or a file of taxonomy IDs to 6 level rank or viceversa "
    )
    print("\n --- When you want to close ExTaxsI just press: " + color.BOLD +
          "CTRL + C" + color.END + " ---\n")

    while True:
        try:
            module = int(input(">> Enter the number of the chosen module: "))
            if module not in (1, 2, 3, 4, 5):
                print("Error, enter a valid choice!")
                continue
            break
        except ValueError:
            print("Wrong input, expected a numeric input, try again.")

    clear()

    if module == 1:
        ETIDB.database_module(None, 0, None, None, [None, None])

    if module == 2:
        ETILSTAT.statistical_module()

    if module == 3:
        taxonomyID_module()

    if module == 4:
        anno = str(input(">>Annotation: "))
        print_data(retrieve_annotation(anno))

    if module == 5:
        print("Oh! You found me out! Here take my money!")
        time.sleep(5)

    print("See you again!!")

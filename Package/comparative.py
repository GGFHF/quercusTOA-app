#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=broad-except
# pylint: disable=invalid-name
# pylint: disable=line-too-long
# pylint: disable=multiple-statements
# pylint: disable=too-many-lines
# pylint: disable=unnecessary-pass

#-------------------------------------------------------------------------------

'''
This file contains the classes related to the comparative genomics of quercusTOA
(Quercus Taxonomy-oriented Annotation).

This software has been developed by:

    GI en Desarrollo de Especies y Comunidades Leñosas (WooSp)
    Dpto. Sistemas y Recursos Naturales
    ETSI Montes, Forestal y del Medio Natural
    Universidad Politecnica de Madrid
    https://github.com/ggfhf/

Licence: GNU General Public Licence Version 3.
'''

#-------------------------------------------------------------------------------

import gzip
import os
import re
import shutil
import subprocess
import sys

from collections import defaultdict
from Bio import Phylo

import matplotlib
import matplotlib.pyplot as plt
from pymsaviz import MsaViz

from PyQt5.QtCore import Qt                      # pylint: disable=no-name-in-module
from PyQt5.QtGui import QCursor                  # pylint: disable=no-name-in-module
from PyQt5.QtGui import QFontMetrics             # pylint: disable=no-name-in-module
from PyQt5.QtGui import QGuiApplication          # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QAbstractItemView    # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QApplication         # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QComboBox            # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QFileDialog          # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QGridLayout          # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QGroupBox            # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QHeaderView          # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QLabel               # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QLineEdit            # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QPlainTextEdit       # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QMessageBox          # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QPushButton          # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QRadioButton         # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QTableWidget         # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QTableWidgetItem     # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QVBoxLayout          # pylint: disable=no-name-in-module
from PyQt5.QtWidgets import QWidget              # pylint: disable=no-name-in-module

import dialogs
import genlib
import sqllib

matplotlib.use('Agg')

#-------------------------------------------------------------------------------

class FormSearchSeqsHomology(QWidget):
    '''
    Class used to search homology of sequences.
    '''

    #---------------

    def __init__(self, parent):
        '''
        Create a class instance.
        '''

        # save parameters in instance variables
        self.parent = parent

        # call the init method of the parent class
        super().__init__()

        # set the dimensions window
        self.window_height = self.parent.WINDOW_HEIGHT - 100
        self.window_width = self.parent.WINDOW_WIDTH - 50

        # set the head and title
        self.head = 'Searh homology of sequences'
        self.title = f'{genlib.get_app_short_name()} - {self.head}'

        # get the dictionary of application configuration
        self.app_config_dict = genlib.get_config_dict(genlib.get_app_config_file())

        # initialize
        self.fasta_source = ''

        # build the graphic user interface of the window
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

        # show the window
        self.show()

    #---------------

    def build_gui(self):
        '''
        Build the graphic user interface of the window.
        '''

        # set the width and height of the window
        self.setFixedSize(self.window_width, self.window_height)

        # move the window at center
        rectangle = self.frameGeometry()
        central_point = QGuiApplication.primaryScreen().availableGeometry().center()
        rectangle.moveCenter(central_point)
        self.move(rectangle.topLeft())

        # get font metrics information
        fontmetrics = QFontMetrics(QApplication.font())

        # create and configure "label_head"
        label_head = QLabel(self.head, alignment=Qt.AlignCenter)
        label_head.setStyleSheet('font: bold 14px; color: black; background-color: lightGray; max-height: 30px')

        # create and configure "label_threads"
        label_threads = QLabel()
        label_threads.setText('Threads')
        label_threads.setFixedWidth(fontmetrics.width('9'*14))

        # create and configure "lineedit_threads"
        self.lineedit_threads  = QLineEdit()
        self.lineedit_threads.setFixedWidth(fontmetrics.width('9'*6))
        self.lineedit_threads.editingFinished.connect(self.check_inputs)

        # create and configure "label_fasta_type"
        label_fasta_type = QLabel()
        label_fasta_type.setText('FASTA type')
        label_fasta_type.setFixedWidth(fontmetrics.width('9'*10))

        # create and configure "combobox_fasta_type"
        self.combobox_fasta_type = QComboBox()
        self.combobox_fasta_type.currentIndexChanged.connect(self.check_inputs)
        self.combobox_fasta_type.setFixedWidth(fontmetrics.width('9'*22))

        # create and configure "radiobutton_seqs_list"
        self.radiobutton_seqs_list = QRadioButton()
        self.radiobutton_seqs_list.setText('FASTA\nsequences')
        self.radiobutton_seqs_list.setFixedWidth(fontmetrics.width('9'*14))
        self.radiobutton_seqs_list.toggled.connect(self.update_radiobutton_state)

        # create and configure "plaintextedit_seqs_list"
        self.plaintextedit_seqs_list = QPlainTextEdit()
        self.plaintextedit_seqs_list.setFixedWidth(fontmetrics.width('9'*75))
        self.plaintextedit_seqs_list.setFixedHeight(fontmetrics.lineSpacing() * 10 + 12)
        self.plaintextedit_seqs_list.textChanged.connect(self.check_inputs)

        # create and configure "radiobutton_fasta_file"
        self.radiobutton_fasta_file = QRadioButton()
        self.radiobutton_fasta_file.setText('FASTA file')
        self.radiobutton_fasta_file.setFixedWidth(fontmetrics.width('9'*14))
        self.radiobutton_fasta_file.toggled.connect(self.update_radiobutton_state)

        # create and configure "lineedit_fasta_file"
        self.lineedit_fasta_file = QLineEdit()
        self.lineedit_fasta_file.editingFinished.connect(self.check_inputs)

        # create and configure "pushbutton_search_fasta_file"
        self.pushbutton_search_fasta_file = QPushButton('Search ...')
        self.pushbutton_search_fasta_file.setToolTip('Search and select the FASTA file.')
        self.pushbutton_search_fasta_file.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_search_fasta_file.clicked.connect(self.pushbutton_search_fasta_file_clicked)

        # create and configure "empty"
        label_empty = QLabel()
        label_empty.setFixedWidth(fontmetrics.width('9'*3))

        # create and configure "gridlayout_data"
        gridlayout_data = QGridLayout()
        gridlayout_data.setRowMinimumHeight(0, 40)
        gridlayout_data.setRowMinimumHeight(1, 40)
        gridlayout_data.setRowMinimumHeight(2, 40)
        gridlayout_data.setColumnStretch(0,1)
        gridlayout_data.setColumnStretch(1,1)
        gridlayout_data.setColumnStretch(2,1)
        gridlayout_data.setColumnStretch(3,1)
        gridlayout_data.setColumnStretch(4,15)
        gridlayout_data.setColumnStretch(5,1)
        gridlayout_data.addWidget(label_threads, 0, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_threads, 0, 1, 1, 1, alignment=Qt.AlignLeft)
        gridlayout_data.addWidget(label_empty, 0, 2, 1, 1)
        gridlayout_data.addWidget(label_fasta_type, 0, 3, 1, 1)
        gridlayout_data.addWidget(self.combobox_fasta_type, 0, 4, 1, 1, alignment=Qt.AlignLeft)
        gridlayout_data.addWidget(self.radiobutton_seqs_list, 1, 0, 1, 1)
        gridlayout_data.addWidget(self.plaintextedit_seqs_list, 1, 1, 1, 5, alignment=Qt.AlignLeft)
        gridlayout_data.addWidget(self.radiobutton_fasta_file, 2, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_fasta_file, 2, 1, 1, 4)
        gridlayout_data.addWidget(self.pushbutton_search_fasta_file, 2, 5, 1, 1)

        # create and configure "groupbox_data"
        groupbox_data = QGroupBox()
        groupbox_data.setObjectName('groupbox_data')
        groupbox_data.setStyleSheet('QGroupBox#groupbox_data {border: 0px;}')
        groupbox_data.setLayout(gridlayout_data)

        # create and configure "pushbutton_execute"
        self.pushbutton_execute = QPushButton('Execute')
        self.pushbutton_execute.setToolTip('Execute the sequences homology search.')
        self.pushbutton_execute.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_execute.clicked.connect(self.pushbutton_execute_clicked)

        # create and configure "pushbutton_close"
        pushbutton_close = QPushButton('Close')
        pushbutton_close.setToolTip('Cancel the run of the sequences homology search and close the window.')
        pushbutton_close.setCursor(QCursor(Qt.PointingHandCursor))
        pushbutton_close.clicked.connect(self.pushbutton_close_clicked)

        # create and configure "gridlayout_buttons"
        gridlayout_buttons = QGridLayout()
        gridlayout_buttons.setColumnStretch(0, 10)
        gridlayout_buttons.setColumnStretch(1, 1)
        gridlayout_buttons.setColumnStretch(2, 1)
        gridlayout_buttons.addWidget(self.pushbutton_execute, 0, 1, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(pushbutton_close, 0, 2, alignment=Qt.AlignCenter)

        # create and configure "groupbox_buttons"
        groupbox_buttons = QGroupBox()
        groupbox_buttons.setObjectName('groupbox_buttons')
        groupbox_buttons.setStyleSheet('QGroupBox#groupbox_buttons {border: 0px;}')
        groupbox_buttons.setLayout(gridlayout_buttons)

        # create and configure "gridlayout_central"
        gridlayout_central = QGridLayout()
        gridlayout_central.setRowStretch(0, 1)
        gridlayout_central.setRowStretch(1, 1)
        gridlayout_central.setRowStretch(2, 1)
        gridlayout_central.setRowStretch(3, 1)
        gridlayout_central.setRowStretch(4, 1)
        gridlayout_central.setColumnStretch(0, 0)
        gridlayout_central.setColumnStretch(1, 1)
        gridlayout_central.setColumnStretch(2, 0)
        gridlayout_central.addWidget(label_head, 0, 1)
        gridlayout_central.addWidget(QLabel(), 1, 1)
        gridlayout_central.addWidget(groupbox_data, 2, 1)
        gridlayout_central.addWidget(QLabel(), 3, 1)
        gridlayout_central.addWidget(groupbox_buttons, 4, 1)

        # create and configure "groupbox_central"
        groupbox_central = QGroupBox()
        groupbox_central.setLayout(gridlayout_central)

        # create and configure "vboxlayout"
        vboxlayout = QVBoxLayout(self)
        vboxlayout.addWidget(groupbox_central)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # set initial value in "lineedit_threads"
        self.lineedit_threads.setText('4')

        # populate data in "combobox_fasta_type"
        self.combobox_fasta_type_populate()

        # set initial value in "plaintextedit_seqs_list"
        self.plaintextedit_seqs_list.setPlaceholderText('')

        # set initial value in "lineedit_fasta_file"
        self.lineedit_fasta_file.setText('')

        # check "radiobutton_seqs_list"
        self.radiobutton_seqs_list.setChecked(True)
        self.fasta_source = 'SEQUENCES'

    #---------------

    def check_inputs(self):
        '''
        Check the content of each input and do the actions linked to its value.
        '''

        # initialize the control variable
        OK = True

        # check "lineedit_threads" when the editing finished
        if not self.lineedit_threads_editing_finished():
            OK = False

        # check "plaintextedit_seqs_list" or "lineedit_fasta_file" when the editing finished
        if not self.plaintextedit_seqs_list_editing_finished() and not self.lineedit_fasta_file_editing_finished():
            OK = False

        # check all inputs are OK
        if OK:
            self.parent.statusBar().showMessage('')
        else:
            self.parent.statusBar().showMessage('WARNING: One or more input values are wrong or empty.')

        # enable "pushbutton_execute"
        if OK and self.lineedit_threads.text() != '' and (self.lineedit_fasta_file.text() != '' or  self.plaintextedit_seqs_list.toPlainText() != ''):
            self.pushbutton_execute.setEnabled(True)
        else:
            self.pushbutton_execute.setEnabled(False)

        # return the control variable
        return OK

    #---------------

    def update_radiobutton_state(self):
        '''
        Check the state of radio buttons and its associated widgets.
        '''

        if self.radiobutton_seqs_list.isChecked():
            self.plaintextedit_seqs_list.setEnabled(True)
            self.lineedit_fasta_file.setText('')
            self.lineedit_fasta_file.setEnabled(False)
            self.pushbutton_search_fasta_file.setEnabled(False)
            self.fasta_source = 'SEQUENCES'

        elif self.radiobutton_fasta_file.isChecked():
            self.plaintextedit_seqs_list.clear()
            self.plaintextedit_seqs_list.setEnabled(False)
            self.lineedit_fasta_file.setEnabled(True)
            self.pushbutton_search_fasta_file.setEnabled(True)
            self.fasta_source = 'FILE'

    #---------------

    def lineedit_threads_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_threads"
        '''

        # initialize the control variable
        OK = True

        # chek if "lineedit_threads" is empty
        if self.lineedit_threads.text() == '':
            OK = False
            self.lineedit_threads.setStyleSheet('background-color: white')

        # chek if "lineedit_threads" is an integer number between 1 and the CPUs number in the computer
        elif self.lineedit_threads.text() != '' and not genlib.check_int(self.lineedit_threads.text(), minimum=1, maximum=os.cpu_count()):
            OK = False
            self.lineedit_threads.setStyleSheet('background-color: red')
            text = f'The value of threads number has to be an integer number between 1 and {os.cpu_count()} (threads available in the computer).'
            QMessageBox.critical(self, self.title, text, buttons=QMessageBox.Ok)

        else:
            self.lineedit_threads.setStyleSheet('background-color: white')

        # return the control variable
        return OK

    #---------------

    def combobox_fasta_type_populate(self):
        '''
        Populate data in "combobox_fasta_type".
        '''

        # populate data in "combobox_fasta_type"
        self.combobox_fasta_type.addItems([genlib.get_fasta_type_transcripts(), genlib.get_fasta_type_proteins()])

        # simultate "combobox_fasta_type" index has changed
        self.combobox_fasta_type_currentIndexChanged()

    #---------------

    def combobox_fasta_type_currentIndexChanged(self):
        '''
        Process the event when an item of "combobox_fasta_type" has been selected.
        '''

        # check the content of inputs
        self.check_inputs()

    #---------------

    def lineedit_fasta_file_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_fasta_file"
        '''

        # initialize the control variable
        OK = True

        # chek if "lineedit_fasta_file" is empty
        if self.lineedit_fasta_file.text() == '':
            self.lineedit_fasta_file.setStyleSheet('background-color: white')
            OK = False

        # chek if "lineedit_fasta_file" exists
        elif not os.path.isfile(self.lineedit_fasta_file.text()):
            self.lineedit_fasta_file.setStyleSheet('background-color: red')
            OK = False

        else:
            self.lineedit_fasta_file.setStyleSheet('background-color: white')

        # return the control variable
        return OK

    #---------------

    def pushbutton_search_fasta_file_clicked(self):
        '''
        Search and select the FASTA file.
        '''

        # search the FASTA file
        (fasta_file, _) = QFileDialog.getOpenFileName(self, f'{self.head} - Selection of FASTA file', os.path.expanduser('~'), "FASTA files (*.fasta *.FASTA *.fas *.FAS *.fa *.FA *.fsa *.FSA *.fna *.FNA *.faa *.FAA);;all (*.*)")

        # set "lineedit_fasta_file" with the FASTA file selected
        if fasta_file != '':
            self.lineedit_fasta_file.setText(fasta_file)

        # check the content of inputs
        self.check_inputs()

    #---------------

    def plaintextedit_seqs_list_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "plaintextedit_seqs_list"
        '''

        # initialize the control variable
        OK = True

        # chek if "plaintextedit_seqs_list" is empty
        if self.plaintextedit_seqs_list.toPlainText() == '':
            OK = False
            self.plaintextedit_seqs_list.setStyleSheet('background-color: white')

        # return the control variable
        return OK

    #---------------

    def pushbutton_execute_clicked(self):
        '''
        Execute the process.
        '''

        # initialize the control variable
        OK = True

        # confirm the process is executed
        text = 'The search of sequences homology is going to be run.\n\nAre you sure to continue?'
        botton = QMessageBox.question(self, self.title, text, buttons=QMessageBox.Yes|QMessageBox.No, defaultButton=QMessageBox.No)
        if botton == QMessageBox.No:
            OK = False

        # execute the process
        if OK:

            # get the threads number
            threads = self.lineedit_threads.text()

            # get the FASTA type
            fasta_type = self.combobox_fasta_type.currentText()

            # get the FASTA file
            fasta_sequences = ''
            fasta_file = ''
            if self.fasta_source == 'SEQUENCES':
                fasta_sequences = self.plaintextedit_seqs_list.toPlainText()
            elif self.fasta_source == 'FILE':
                fasta_file = self.lineedit_fasta_file.text()
                if sys.platform.startswith('win32'):
                    fasta_file = genlib.windows_path_2_wsl_path(fasta_file)

            # create and execute "DialogProcess"
            process = dialogs.DialogProcess(self, self.head, self.search_seqs_homology, threads, fasta_type, self.fasta_source, fasta_sequences, fasta_file)
            process.exec()

        # close the windows
        # -- if OK:
        # --     self.pushbutton_close_clicked()

    #---------------

    def pushbutton_close_clicked(self):
        '''
        Close the window.
        '''

        self.parent.current_subwindow = None
        self.close()
        self.parent.set_background_image()

   #---------------

    def search_seqs_homology(self, process, threads, fasta_type, fasta_source, fasta_sequences, fasta_file):
        '''
        Run a search of sequences homology.
        '''

        # initialize the control variable
        OK = True

        # warn that the log window does not have to be closed
        process.write('Do not close this window, please wait!\n')

        # warn that the requirements are being verified
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Checking process requirements ...\n')

        # check the quercusTOA config file
        if OK:
            if not os.path.isfile(genlib.get_app_config_file()):
                process.write(f'*** ERROR: The {genlib.get_app_short_name()} config file does not exist. Please, recreate it.\n')
                OK = False

        # warn that the requirements are OK
        if OK:
            process.write('Process requirements are OK.\n')

        # determine the temporal directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Determining the temporal directory ...\n')
            temp_dir = genlib.get_temp_dir()
            command = f'mkdir -p {temp_dir}'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write(f'The directory path is {temp_dir}.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # determine the run directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Determining the run directory ...\n')
            result_dir = self.app_config_dict['Environment parameters']['result_dir']
            current_run_dir = genlib.get_current_run_dir(result_dir, genlib.get_result_run_subdir(), genlib.get_process_search_seqs_homology_code())
            command = f'mkdir -p {current_run_dir}'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write(f'The directory path is {current_run_dir}.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # build the script
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            script_name = f'{genlib.get_process_search_seqs_homology_code()}-process.sh'
            process.write(f'Building the process script {script_name} ...\n')
            (OK, _) = self.build_search_seqs_homology_script(temp_dir, script_name, current_run_dir, threads, fasta_type, fasta_source, fasta_sequences, fasta_file)
            if OK:
                process.write('The file is built.\n')
            else:
                process.write('*** ERROR: The file could not be built.\n')

        # copy the script to the current run directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Copying the script {script_name} to the directory {current_run_dir} ...\n')
            command = f'cp {temp_dir}/{script_name} {current_run_dir}; [ $? -eq 0 ] &&  exit 0 || exit 1'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write('The file is copied.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')

        # set run permision to the script
        if OK and not sys.platform.startswith('win32'):
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Setting on the run permision of {script_name} ...\n')
            command = f'chmod u+x {current_run_dir}/{script_name}'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write('The run permision is set.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # build the starter script in the temporal directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            starter_name = f'{genlib.get_process_search_seqs_homology_code()}-process-starter.sh'
            process.write(f'Building the starter script {starter_name} ...\n')
            (OK, _) = genlib.build_starter(temp_dir, starter_name, script_name, current_run_dir)
            if OK:
                process.write('The file is built.\n')
            if not OK:
                process.write('***ERROR: The file could not be built.\n')

        # copy the starter script to the current run directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Copying the starter {starter_name} to the directory {current_run_dir} ...\n')
            command = f'cp {temp_dir}/{starter_name} {current_run_dir}; [ $? -eq 0 ] &&  exit 0 || exit 1'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write('The file is copied.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')

        # set run permision to the starter script
        if OK and not sys.platform.startswith('win32'):
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Setting on the run permision of {starter_name} ...\n')
            command = f'chmod u+x {current_run_dir}/{starter_name}'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write('The run permision is set.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # submit the starter
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Submitting the starter {starter_name} ...\n')
            command = f'{current_run_dir}/{starter_name} &'
            rc = genlib.run_command(command, process, is_script=True)
            if rc == 0:
                process.write('The script is submitted.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # warn that the log window can be closed
        process.write(f'{genlib.get_separator()}\n')
        process.write('You can close this window now.\n')

        # return the control variable
        return OK

    #---------------

    def build_search_seqs_homology_script(self, directory, script_name, current_run_dir, threads, fasta_type, fasta_source, fasta_sequences, fasta_file):
        '''
        Build the script to run a homology relationships pipeline.
        '''

        # initialize the control variable and error list
        OK = True
        error_list = []

        # get the Miniforge3 directory and its bin subdirectory
        miniforge3_dir = ''
        if sys.platform.startswith('win32'):
            miniforge3_dir = genlib.get_miniforge3_dir_in_wsl()
        elif sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
            miniforge3_dir = genlib.get_miniforge3_current_dir()
        miniforge3_bin_dir = f'{miniforge3_dir}/bin'

        # get items from dictionary of application configuration
        app_dir = self.app_config_dict['Environment parameters']['app_dir']
        sequences_db_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['sequences_db_path']
        functional_annotations_db_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['functional_annotations_db_path']
        comparative_genomics_db_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['comparative_genomics_db_path']
        quercus_sequence_file = self.app_config_dict[f'{genlib.get_app_short_name()} database']['quercus_sequence_file']
        quercus_blastplus_db_name = self.app_config_dict[f'{genlib.get_app_short_name()} database']['quercus_blastplus_db_name']
        quercus_blastplus_db_dir = self.app_config_dict[f'{genlib.get_app_short_name()} database']['quercus_blastplus_db_dir']
        quercus_diamond_db_name = self.app_config_dict[f'{genlib.get_app_short_name()} database']['quercus_diamond_db_name']
        quercus_diamond_db_dir = self.app_config_dict[f'{genlib.get_app_short_name()} database']['quercus_diamond_db_dir']

        # set the CodAn output directory
        codan_output_dir = f'{current_run_dir}/codan_output'

        # set the parameters file
        params_file = f'{current_run_dir}/{genlib.get_params_file_name()}'

        # set the temporal directories
        temp_dir = f'{current_run_dir}/temp'
        temp_liftoff_dir = f'{temp_dir}/temp-liftoff'

        # when the FASTA source is sequences set the FASTA file
        if fasta_source == 'SEQUENCES':
            fasta_file = f'{temp_dir}/analysis-sequences.fasta'

        # set the CSV file with alignments
        blastp_alignment_file = f'{temp_dir}/{genlib.get_blastp_clade_alignment_file_name()}'

        # set the homology relationships file
        homology_relationships_file = f'{current_run_dir}/{genlib.get_homology_relationships_file_name()}'

        # set outdir directory of the protein FASTA files corresponding to homology relationships and their alignments
        seqs_alignments_dir = f'{current_run_dir}/seqs-alignments'

        # set the default parameters of the alignment tool
        codan_model = self.app_config_dict['CodAn models']['codan_full_plants_model_dir']
        alignment_tool = genlib.get_diamond_name()
        evalue = 1E-6
        max_target_seqs = 1
        max_hsps = 1
        qcov_hsp_perc = 0.0
        other_parameters = 'NONE'

        # set the script path
        script_path = f'{directory}/{script_name}'

        # write the script
        try:
            with open(script_path, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
                file_id.write( '#!/bin/bash\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write(f'export PATH={miniforge3_bin_dir}:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin\n')
                file_id.write( 'SEP="#########################################"\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write(f'STATUS_DIR={genlib.get_status_dir(current_run_dir)}\n')
                file_id.write(f'SCRIPT_STATUS_OK={genlib.get_status_ok(current_run_dir)}\n')
                file_id.write(f'SCRIPT_STATUS_WRONG={genlib.get_status_wrong(current_run_dir)}\n')
                file_id.write( 'mkdir -p $STATUS_DIR\n')
                file_id.write(f'mkdir -p {temp_dir}\n')
                file_id.write(f'mkdir -p {temp_liftoff_dir}\n')
                file_id.write(f'mkdir -p {seqs_alignments_dir}\n')
                file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function init\n')
                file_id.write( '{\n')
                file_id.write( '    INIT_DATETIME=`date +%s`\n')
                file_id.write( '    FORMATTED_INIT_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME."\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write('function save_params\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Saving parameters ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/save-parms.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                file_id.write(f'        echo "[Search parameters]" > {params_file}\n')
                file_id.write(f'        echo "fasta_type = {fasta_type}" >> {params_file}\n')
                file_id.write(f'        echo "fasta_file = {fasta_file}" >> {params_file}\n')
                file_id.write(f'        echo "codan_model = {codan_model}" >> {params_file}\n')
                file_id.write(f'        echo "alignment_tool = {alignment_tool}" >> {params_file}\n')
                file_id.write(f'        echo "evalue = {evalue}" >> {params_file}\n')
                file_id.write(f'        echo "max_target_seqs = {max_target_seqs}" >> {params_file}\n')
                file_id.write(f'        echo "max_hsps = {max_hsps}" >> {params_file}\n')
                file_id.write(f'        echo "qcov_hsp_perc = {qcov_hsp_perc}" >> {params_file}\n')
                file_id.write(f'        echo "other_parameters = {other_parameters}" >> {params_file}\n')
                file_id.write( '        RC=$?\n')
                file_id.write( '        if [ $RC -ne 0 ]; then manage_error echo $RC; fi\n')
                file_id.write( '        echo "Parameters are saved."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                if fasta_source == 'SEQUENCES':
                    file_id.write( '#-------------------------------------------------------------------------------\n')
                    file_id.write('function save_fasta_squences\n')
                    file_id.write( '{\n')
                    file_id.write( '    echo "$SEP"\n')
                    file_id.write( '    echo "Saving fasta sequences ..."\n')
                    file_id.write(f'    cd {current_run_dir}\n')
                    file_id.write( '    STEP_STATUS=$STATUS_DIR/save-fasta-sequences.ok\n')
                    file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                    file_id.write( '        echo "This step was previously run."\n')
                    file_id.write( '    else\n')
                    file_id.write(f'        echo "{fasta_sequences}" > {fasta_file}\n')
                    file_id.write( '        RC=$?\n')
                    file_id.write( '        if [ $RC -ne 0 ]; then manage_error echo $RC; fi\n')
                    file_id.write( '        echo "Sequences are saved."\n')
                    file_id.write( '        touch $STEP_STATUS\n')
                    file_id.write( '    fi\n')
                    file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function predict_orfs\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Predicting ORFs and getting peptide sequences ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/predict-orfs.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                if fasta_type ==  genlib.get_fasta_type_transcripts():
                    file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_codan_environment()}\n')
                    file_id.write( '        MODELS_DIR=`echo $CONDA_PREFIX`/models\n')
                    file_id.write( '        /usr/bin/time \\\n')
                    file_id.write( '            codan.py \\\n')
                    file_id.write(f'                --cpu={threads} \\\n')
                    file_id.write(f'                --model=$MODELS_DIR/{codan_model} \\\n')
                    file_id.write(f'                --transcripts={fasta_file} \\\n')
                    file_id.write(f'                --output={codan_output_dir}\n')
                    file_id.write( '        RC=$?\n')
                    file_id.write( '        if [ $RC -ne 0 ]; then manage_error codan.py $RC; fi\n')
                    if codan_model.endswith('_partial'):
                        file_id.write( '        /usr/bin/time \\\n')
                        file_id.write( '            TranslatePartial.py \\\n')
                        file_id.write(f'                {codan_output_dir}/ORF_sequences.fasta \\\n')
                        file_id.write(f'                {codan_output_dir}/PEP_sequences.fa\n')
                        file_id.write( '        RC=$?\n')
                        file_id.write( '        if [ $RC -ne 0 ]; then manage_error TranslatePartial.py $RC; fi\n')
                    file_id.write( '        conda deactivate\n')
                    file_id.write( '        echo "ORFs are predicted and peptide sequences are gotten."\n')
                elif fasta_type ==  genlib.get_fasta_type_proteins():
                    file_id.write( '        echo "This step is not run with a proteins file."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function align_peptides_2_alignment_tool_quercus_db\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write(f'    echo "Aligning peptides to the {alignment_tool} Quercus database ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/align-peptides-2-alignment-tool-quercus-db.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                if fasta_type ==  genlib.get_fasta_type_transcripts():
                    file_id.write(f'        PEPTIDE_FILE={codan_output_dir}/PEP_sequences.fa\n')
                elif fasta_type ==  genlib.get_fasta_type_proteins():
                    file_id.write(f'        PEPTIDE_FILE={fasta_file}\n')
                if alignment_tool == genlib.get_blastplus_name():
                    file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_blastplus_environment()}\n')
                    file_id.write(f'        export BLASTDB={quercus_blastplus_db_dir}\n')
                    file_id.write( '        /usr/bin/time \\\n')
                    file_id.write( '            blastp \\\n')
                    file_id.write(f'                -num_threads {threads} \\\n')
                    file_id.write(f'                -db {quercus_blastplus_db_name} \\\n')
                    file_id.write( '                -query $PEPTIDE_FILE \\\n')
                    file_id.write(f'                -evalue {evalue} \\\n')
                    file_id.write(f'                -max_target_seqs {max_target_seqs} \\\n')
                    file_id.write(f'                -max_hsps {max_hsps} \\\n')
                    file_id.write(f'                -qcov_hsp_perc {qcov_hsp_perc} \\\n')
                    file_id.write( '                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\\n')
                    if other_parameters.upper() != 'NONE':
                        parameter_list = [x.strip() for x in other_parameters.split(';')]
                        for parameter in parameter_list:
                            if parameter.find('=') > 0:
                                pattern = r'^--(.+)=(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                parameter_value = mo.group(2).strip()
                                file_id.write(f'                -{parameter_name} {parameter_value} \\\n')
                            else:
                                pattern = r'^--(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                file_id.write(f'                -{parameter_name} \\\n')
                    file_id.write(f'                -out {blastp_alignment_file}\n')
                    file_id.write( '        RC=$?\n')
                    file_id.write( '        if [ $RC -ne 0 ]; then manage_error blastp $RC; fi\n')
                    file_id.write( '        conda deactivate\n')
                elif alignment_tool == genlib.get_diamond_name():
                    file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_diamond_environment()}\n')
                    file_id.write( '        /usr/bin/time \\\n')
                    file_id.write( '            diamond blastp \\\n')
                    file_id.write(f'                --threads {threads} \\\n')
                    file_id.write(f'                --db {quercus_diamond_db_dir}/{quercus_diamond_db_name} \\\n')
                    file_id.write( '                --query $PEPTIDE_FILE \\\n')
                    file_id.write(f'                --evalue {evalue} \\\n')
                    file_id.write(f'                --max-target-seqs {max_target_seqs} \\\n')
                    file_id.write(f'                --max-hsps {max_hsps} \\\n')
                    file_id.write( '                --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\\n')
                    if other_parameters.upper() != 'NONE':
                        parameter_list = [x.strip() for x in other_parameters.split(';')]
                        for parameter in parameter_list:
                            if parameter.find('=') > 0:
                                pattern = r'^--(.+)=(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                parameter_value = mo.group(2).strip()
                                file_id.write(f'                --{parameter_name} {parameter_value} \\\n')
                            else:
                                pattern = r'^--(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                file_id.write(f'                --{parameter_name} \\\n')
                    file_id.write(f'                --out {blastp_alignment_file}\n')
                    file_id.write( '        RC=$?\n')
                    file_id.write( '        if [ $RC -ne 0 ]; then manage_error diamond-blastp $RC; fi\n')
                    file_id.write( '        conda deactivate\n')
                file_id.write( '        echo "Alignment is done."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write('function get_cluster_homology_relationships\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Getting cluster homology relationships ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/get-cluster-homology-relationships.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_quercustoa_env_code()}\n')
                file_id.write( '        /usr/bin/time \\\n')
                file_id.write(f'            {app_dir}/get-cluster-homology-relationships.py \\\n')
                file_id.write(f'                --comparative-db={comparative_genomics_db_path} \\\n')
                file_id.write(f'                --annotations-db={functional_annotations_db_path} \\\n')
                file_id.write(f'                --blastp-alignments={blastp_alignment_file} \\\n')
                file_id.write(f'                --homology={homology_relationships_file} \\\n')
                file_id.write( '                --verbose=N \\\n')
                file_id.write( '                --trace=N\n')
                file_id.write( '        RC=$?\n')
                file_id.write( '        if [ $RC -ne 0 ]; then manage_error get-cluster-homology-relationships.py $RC; fi\n')
                file_id.write( '        conda deactivate\n')
                file_id.write( '        echo "Cluster homology relationships are gotten."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write('function get_protein_fasta_files\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Getting protein FASTA files ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/get-protein-fasta-files.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_quercustoa_env_code()}\n')
                file_id.write( '        /usr/bin/time \\\n')
                file_id.write(f'            {app_dir}/get-protein-fasta-files.py \\\n')
                file_id.write(f'                --sequences-db={sequences_db_path} \\\n')
                file_id.write(f'                --homology={homology_relationships_file} \\\n')
                file_id.write(f'                --blastp-alignments={blastp_alignment_file} \\\n')
                if fasta_type ==  genlib.get_fasta_type_transcripts():
                    file_id.write( '                --analysis=NONE \\\n')
                elif fasta_type ==  genlib.get_fasta_type_proteins():
                    file_id.write(f'                --analysis={fasta_file} \\\n')
                file_id.write(f'                --consensus={quercus_sequence_file} \\\n')
                file_id.write(f'                --outdir={seqs_alignments_dir} \\\n')
                file_id.write( '                --verbose=N \\\n')
                file_id.write( '                --trace=N\n')
                file_id.write( '        RC=$?\n')
                file_id.write( '        if [ $RC -ne 0 ]; then manage_error get-cluster-homology-relationships.py $RC; fi\n')
                file_id.write( '        conda deactivate\n')
                file_id.write( '        echo "Protein FASTA files are gotten."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write('function get_alignment_files\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Getting alignment files and their plots ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/get_alignment_files.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_mafft_environment()}\n')
                file_id.write(f'        FASTA_FILE_LIST={seqs_alignments_dir}/fasta_file_list.txt\n')
                file_id.write(f'        ls {seqs_alignments_dir}/*.fasta > $FASTA_FILE_LIST\n')
                file_id.write( '        while read FASTA_FILE; do\n')
                file_id.write( '            ALIGNMENT_FILE=`echo $FASTA_FILE | sed "s/.fasta/.aln/g"`\n')
                file_id.write( '            if [[ "$FASTA_FILE" == *homologous-proteins.fasta ]]; then\n')
                file_id.write( '                TREE=Y\n')
                file_id.write( '            else\n')
                file_id.write( '                TREE=N\n')
                file_id.write( '            fi\n')
                file_id.write( '            /usr/bin/time \\\n')
                file_id.write(f'                {app_dir}/align-fasta-seqs.py \\\n')
                file_id.write( '                    --seqs=$FASTA_FILE \\\n')
                file_id.write( '                    --tree=$TREE \\\n')
                file_id.write( '                    --verbose=N \\\n')
                file_id.write( '                    --trace=N\n')
                file_id.write( '            RC=$?\n')
                file_id.write( '            if [ $RC -ne 0 ]; then manage_error align-fasta-seqs.py $RC; fi\n')
                file_id.write( '        done < $FASTA_FILE_LIST\n')
                file_id.write( '        conda deactivate\n')
                file_id.write( '        echo "FASTA files are aligned."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function end\n')
                file_id.write( '{\n')
                file_id.write( '    END_DATETIME=`date +%s`\n')
                file_id.write( '    FORMATTED_END_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`\n')
                file_id.write( '    calculate_duration\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Script ended OK at $FORMATTED_END_DATETIME with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    touch $SCRIPT_STATUS_OK\n')
                file_id.write( '    exit 0\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function manage_error\n')
                file_id.write( '{\n')
                file_id.write( '    END_DATETIME=`date +%s`\n')
                file_id.write( '    FORMATTED_END_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`\n')
                file_id.write( '    calculate_duration\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "ERROR: $1 returned error $2"\n')
                file_id.write( '    echo "Script ended WRONG at $FORMATTED_END_DATETIME with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    touch $SCRIPT_STATUS_WRONG\n')
                file_id.write( '    exit 3\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function calculate_duration\n')
                file_id.write( '{\n')
                file_id.write( '    DURATION=`expr $END_DATETIME - $INIT_DATETIME`\n')
                file_id.write( '    HH=`expr $DURATION / 3600`\n')
                file_id.write( '    MM=`expr $DURATION % 3600 / 60`\n')
                file_id.write( '    SS=`expr $DURATION % 60`\n')
                file_id.write( '    FORMATTED_DURATION=`printf "%03d:%02d:%02d\\n" $HH $MM $SS`\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'init\n')
                file_id.write( 'save_params\n')
                if fasta_source == 'SEQUENCES':
                    file_id.write( 'save_fasta_squences\n')
                file_id.write( 'predict_orfs\n')
                file_id.write( 'align_peptides_2_alignment_tool_quercus_db\n')
                file_id.write( 'get_cluster_homology_relationships\n')
                file_id.write( 'get_protein_fasta_files\n')
                file_id.write( 'get_alignment_files\n')
                file_id.write( 'end\n')
        except Exception as e:
            error_list.append(f'*** EXCEPTION: "{e}".')
            error_list.append(f'*** ERROR: The file {script_path} is not created.')
            OK = False

        # return the control variable and error list
        return (OK, error_list)

    #---------------

#-------------------------------------------------------------------------------

class FormRestartHomologySearch(QWidget):
    '''
    Class used to restart a sequences homology search.
    '''

    #---------------

    def __init__(self, parent):
        '''
        Create a class instance.
        '''

        # save parameters in instance variables
        self.parent = parent

        # call the init method of the parent class
        super().__init__()

        # set the dimensions window
        self.window_height = self.parent.WINDOW_HEIGHT - 100
        self.window_width = self.parent.WINDOW_WIDTH - 50

        # set the head and title
        self.head = 'Restart a sequences homology search'
        self.title = f'{genlib.get_app_short_name()} - {self.head}'

        # get the dictionary of application configuration
        self.app_config_dict = genlib.get_config_dict(genlib.get_app_config_file())

        # build the graphic user interface of the window
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

        # show the window
        self.show()

    #---------------

    def build_gui(self):
        '''
        Build the graphic user interface of the window.
        '''

        # set the width and height of the window
        self.setFixedSize(self.window_width, self.window_height)

        # move the window at center
        rectangle = self.frameGeometry()
        central_point = QGuiApplication.primaryScreen().availableGeometry().center()
        rectangle.moveCenter(central_point)
        self.move(rectangle.topLeft())

        # get font metrics information
        fontmetrics = QFontMetrics(QApplication.font())

        # create and configure "label_head"
        label_head = QLabel(self.head, alignment=Qt.AlignCenter)
        label_head.setStyleSheet('font: bold 14px; color: black; background-color: lightGray; max-height: 30px')

        # create and configure "tablewidget"
        self.tablewidget = QTableWidget()
        self.tablewidget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tablewidget.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.column_name_list = ['Process', 'Result dataset', 'Date', 'Time', 'Status']
        self.tablewidget.setColumnCount(len(self.column_name_list))
        self.tablewidget.setHorizontalHeaderLabels(self.column_name_list)
        self.tablewidget.setColumnWidth(0, 230)
        self.tablewidget.setColumnWidth(1, 280)
        self.tablewidget.setColumnWidth(2, 85)
        self.tablewidget.setColumnWidth(3, 70)
        self.tablewidget.setColumnWidth(4, 90)
        self.tablewidget.verticalHeader().setVisible(True)
        self.tablewidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.tablewidget.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tablewidget.currentCellChanged.connect(self.tablewidget_currentCellChanged)
        self.tablewidget.cellClicked.connect(self.tablewidget_cellClicked)
        self.tablewidget.cellDoubleClicked.connect(self.tablewidget_cellDoubleClicked)

        # create and configure "label_fasta_type"
        label_fasta_type = QLabel()
        label_fasta_type.setText('FASTA type')
        label_fasta_type.setFixedWidth(fontmetrics.width('9'*10))

        # create and configure "lineedit_fasta_type"
        self.lineedit_fasta_type = QLineEdit()
        self.lineedit_fasta_type.editingFinished.connect(self.check_inputs)
        self.lineedit_fasta_type.setDisabled(True)

        # create and configure "label_fasta_file"
        label_fasta_file = QLabel()
        label_fasta_file.setText('FASTA file')
        label_fasta_file.setFixedWidth(fontmetrics.width('9'*20))

        # create and configure "lineedit_fasta_file"
        self.lineedit_fasta_file = QLineEdit()
        self.lineedit_fasta_file.editingFinished.connect(self.check_inputs)
        self.lineedit_fasta_file.setDisabled(True)

        # create and configure "gridlayout_data"
        gridlayout_data = QGridLayout()
        gridlayout_data.setColumnStretch(0, 1)
        gridlayout_data.setColumnStretch(1, 1)
        gridlayout_data.addWidget(self.tablewidget, 0, 0, 1, 11)
        gridlayout_data.addWidget(label_fasta_type, 1, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_fasta_type, 1, 1, 1, 1)
        gridlayout_data.addWidget(label_fasta_file, 2, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_fasta_file, 2, 1, 1, 10)

        # create and configure "groupbox_data"
        groupbox_data = QGroupBox()
        groupbox_data.setObjectName('groupbox_data')
        groupbox_data.setStyleSheet('QGroupBox#groupbox_data {border: 0px;}')
        groupbox_data.setLayout(gridlayout_data)

        # create and configure "pushbutton_refresh"
        self.pushbutton_refresh = QPushButton('Refresh')
        self.pushbutton_refresh.setToolTip('Update the process list.')
        self.pushbutton_refresh.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_refresh.clicked.connect(self.pushbutton_refresh_clicked)

        # create and configure "pushbutton_execute"
        self.pushbutton_execute = QPushButton('Execute')
        self.pushbutton_execute.setToolTip('Restart the process selected.')
        self.pushbutton_execute.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_execute.clicked.connect(self.pushbutton_execute_clicked)

        # create and configure "pushbutton_close"
        pushbutton_close = QPushButton('Close')
        pushbutton_close.setToolTip('Close the window.')
        pushbutton_close.setCursor(QCursor(Qt.PointingHandCursor))
        pushbutton_close.clicked.connect(self.pushbutton_close_clicked)

        # create and configure "gridlayout_buttons"
        gridlayout_buttons = QGridLayout()
        gridlayout_buttons.setColumnStretch(0, 15)
        gridlayout_buttons.setColumnStretch(1, 1)
        gridlayout_buttons.setColumnStretch(2, 1)
        gridlayout_buttons.setColumnStretch(3, 1)
        gridlayout_buttons.addWidget(self.pushbutton_refresh, 0, 1, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(self.pushbutton_execute, 0, 2, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(pushbutton_close, 0, 3, alignment=Qt.AlignCenter)

        # create and configure "groupbox_buttons"
        groupbox_buttons = QGroupBox()
        groupbox_buttons.setObjectName('groupbox_buttons')
        groupbox_buttons.setStyleSheet('QGroupBox#groupbox_buttons {border: 0px;}')
        groupbox_buttons.setLayout(gridlayout_buttons)

        # create and configure "gridlayout_central"
        gridlayout_central = QGridLayout()
        gridlayout_central.setRowStretch(0, 1)
        gridlayout_central.setRowStretch(1, 1)
        gridlayout_central.setRowStretch(2, 10)
        gridlayout_central.setRowStretch(3, 1)
        gridlayout_central.setColumnStretch(0, 1)
        gridlayout_central.addWidget(label_head, 0, 0)
        gridlayout_central.addWidget(QLabel(), 1, 0)
        gridlayout_central.addWidget(groupbox_data, 2, 0)
        gridlayout_central.addWidget(groupbox_buttons, 3, 0)

        # create and configure "groupbox_central"
        groupbox_central = QGroupBox()
        groupbox_central.setLayout(gridlayout_central)

        # create and configure "vboxlayout"
        vboxlayout = QVBoxLayout(self)
        vboxlayout.addWidget(groupbox_central)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load data in "tablewidget"
        self.load_tablewidget()

        # initialize "lineedit_fasta_type"
        self.lineedit_fasta_type.setText('')

        # initialize "lineedit_fasta_file"
        self.lineedit_fasta_file.setText('')

    #---------------

    def check_inputs(self):
        '''
        Check the content of each input and do the actions linked to its value.
        '''

        # initialize the control variable
        OK = True

        # get the list of rows selected
        row_list = []
        for idx in self.tablewidget.selectionModel().selectedIndexes():
            row_list.append(idx.row())
        row_list = list(set(row_list))

        # check "lineedit_fasta_type" when the editing finished
        if not self.lineedit_fasta_type_editing_finished():
            OK = False

        # check "lineedit_fasta_file" when the editing finished
        if not self.lineedit_fasta_file_editing_finished():
            OK = False

        # check all inputs are OK
        if OK:
            self.parent.statusBar().showMessage('')
        else:
            self.parent.statusBar().showMessage('There are one or more inputs without data or with wrong value.')

        # enable "pushbutton_execute"
        if OK and len(row_list) == 1 and self.lineedit_fasta_type.text() != '' and self.lineedit_fasta_file.text() != '':
            self.pushbutton_execute.setEnabled(True)
        else:
            self.pushbutton_execute.setEnabled(False)
            OK = False

        # return the control variable
        return OK

    #---------------

    def tablewidget_currentCellChanged(self, row, _):
        '''
        Perform necessary actions after changing the current "tablewidget" cell.
        '''

        # check if there is a row selected
        if row >= 0:

            # get the dictionary of application configuration
            result_dir = self.app_config_dict['Environment parameters']['result_dir']

            # get the result dataset corresponding to the row clicked
            result_dataset = self.tablewidget.item(row, 1).text()

            # get the parameters file corresponding to the result dataset
            params_file = f'{result_dir}{os.sep}{genlib.get_result_run_subdir()}{os.sep}{result_dataset}{os.sep}{genlib.get_params_file_name()}'
            if sys.platform.startswith('win32'):
                params_file = genlib.wsl_path_2_windows_path(params_file)

            # get the dictionary of parameters
            params_dict = genlib.get_config_dict(params_file)

            # set the FASTA file in "fasta_type"
            fasta_type = params_dict['Search parameters']['fasta_type']
            self.lineedit_fasta_type.setText(fasta_type)

            # set the FASTA file in "fasta_file"
            fasta_file = params_dict['Search parameters']['fasta_file']
            if sys.platform.startswith('win32'):
                fasta_file = genlib.wsl_path_2_windows_path(fasta_file)
            self.lineedit_fasta_file.setText(fasta_file)

        # check the content of inputs
        self.check_inputs()

    #---------------

    def tablewidget_cellClicked(self, _, __):
        '''
        Perform necessary actions after clicking on a "tablewidget" cell.
        '''

        # check the content of inputs
        self.check_inputs()

    #---------------

    def tablewidget_cellDoubleClicked(self, _, __):
        '''
        Perform necessary actions after double clicking on "tablewidget" cell.
        '''

        # check the content of inputs
        OK = self.check_inputs()

        # if inputs are OK, simulate a click on "pushbutton_execute"
        if OK:
            self.pushbutton_execute_clicked()

    #---------------

    def lineedit_fasta_type_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_fasta_type"
        '''

        # initialize the control variable
        OK = True

        # return the control variable
        return OK

    #---------------

    def lineedit_fasta_file_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_fasta_file"
        '''

        # initialize the control variable
        OK = True

        # return the control variable
        return OK

    #---------------

    def pushbutton_refresh_clicked(self):
        '''
        Refresh "tablewidget" and initialize other inputs.
        '''

        # check the content of inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

    #---------------

    def pushbutton_execute_clicked(self):
        '''
        Execute the process.
        '''

        # initialize the control variable
        OK = True

        # get the list of rows selected
        row_list = []
        for idx in self.tablewidget.selectionModel().selectedIndexes():
            row_list.append(idx.row())
        row_list = list(set(row_list))

        # check if there is only a row selected
        if len(row_list) != 1:
            title = f'{genlib.get_app_short_name()} - {self.head}'
            text = 'One row has to be selected.'
            QMessageBox.critical(self, title, text, buttons=QMessageBox.Ok)
            OK = False

        # confirm the process is executed
        if OK:
            text = 'The sequences homology search is going to be restarted.\n\nAre you sure to continue?'
            botton = QMessageBox.question(self, self.title, text, buttons=QMessageBox.Yes|QMessageBox.No, defaultButton=QMessageBox.No)
            if botton == QMessageBox.No:
                OK = False

        # execute the process
        if OK:

            # get the result dataset
            result_dataset = self.tablewidget.item(row_list[0], 1).text()

            # create and execute "DialogProcess"
            process = dialogs.DialogProcess(self, self.head, self.restart_seqs_homology_search, result_dataset)
            process.exec()

        # close the windows
        # -- if OK:
        # --     self.pushbutton_close_clicked()

    #---------------

    def pushbutton_close_clicked(self):
        '''
        Close the window.
        '''

        self.parent.current_subwindow = None
        self.close()
        self.parent.set_background_image()

    #---------------

    def load_tablewidget(self):
        '''
        Load data in "tablewidget".
        '''

        # get the result directory
        result_dir = self.app_config_dict['Environment parameters']['result_dir']

        # set the type, name and code of the annotation pipeline datasets
        process_type = genlib.get_result_run_subdir()
        process_name = genlib.get_process_search_seqs_homology_name()
        process_code = genlib.get_process_id(process_name)

        # get the process dictionary
        process_dict = genlib.get_process_dict()

        # get the log directory
        log_dir = f'{result_dir}/{process_type}'
        if sys.platform.startswith('win32'):
            log_dir = genlib.wsl_path_2_windows_path(log_dir)

        # set the command to get the result datasets of annotation pipeline in the log directory
        command = ''
        if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
            if process_name == 'all':
                command = f'ls -d {log_dir}/*  | xargs -n 1 basename'
            else:
                command = f'ls -d {log_dir}/{process_code}-*  | xargs -n 1 basename'
        elif sys.platform.startswith('win32'):
            log_dir = log_dir.replace('/', '\\')
            if process_name == 'all':
                command = f'dir /a:d /b {log_dir}\\*'
            else:
                command = f'dir /a:d /b {log_dir}\\{process_code}-*'

        # run the command to get the result datasets of enrichment analysis in the log directory
        output = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=False)

        # initialize the result dataset dictionary
        result_dataset_dict = {}

        # build the result dataset dictionary
        for line in output.stdout.split('\n'):
            if line != '':


                # get data
                result_dataset_id = line.strip()
                try:
                    pattern = r'^(.+)\-(.+)\-(.+)$'
                    mo = re.search(pattern, result_dataset_id)
                    process_code = mo.group(1).strip()
                    process_name = process_dict[process_code]['name']
                    yymmdd = mo.group(2)
                    hhmmss = mo.group(3)
                    date = f'20{yymmdd[:2]}-{yymmdd[2:4]}-{yymmdd[4:]}'
                    time = f'{hhmmss[:2]}:{hhmmss[2:4]}:{hhmmss[4:]}'
                except:    # pylint: disable=bare-except
                    process_name = 'unknown process'
                    date = '0000-00-00'
                    time = '00:00:00'

                # determine the status
                status_ok = os.path.isfile(genlib.get_status_ok(os.path.join(log_dir, result_dataset_id)))
                status_wrong = os.path.isfile(genlib.get_status_wrong(os.path.join(log_dir, result_dataset_id)))
                status = ''
                if status_ok and not status_wrong:
                    status = 'OK'
                elif not status_ok and status_wrong:
                    status = 'wrong'
                elif not status_ok and not status_wrong:
                    status = 'not finished'
                elif status_ok and status_wrong:
                    status = 'undetermined'

                # insert data in the dictionary when the the status is wrong
                if status == 'wrong':
                    key = f'{process_name}-{result_dataset_id}'
                    result_dataset_dict[key] = {'process': process_name, 'result_dataset_id': result_dataset_id, 'date': date, 'time': time, 'status': status}

        # initialize "tablewidget"
        self.tablewidget.clearContents()

        # set the rows number of "tablewidget"
        self.tablewidget.setRowCount(len(result_dataset_dict))

        # load data in "tablewidget" for the OK result datasets of annotation pipeline
        if not result_dataset_dict:
            text = 'There is no run ended wrong.'
            QMessageBox.warning(self, self.title, text, buttons=QMessageBox.Ok)
        else:
            row = 0
            for key in sorted(result_dataset_dict.keys()):
                self.tablewidget.setItem(row, 0, QTableWidgetItem(result_dataset_dict[key]['process']))
                self.tablewidget.setItem(row, 1, QTableWidgetItem(result_dataset_dict[key]['result_dataset_id']))
                self.tablewidget.setItem(row, 2, QTableWidgetItem(result_dataset_dict[key]['date']))
                self.tablewidget.setItem(row, 3, QTableWidgetItem(result_dataset_dict[key]['time']))
                self.tablewidget.setItem(row, 4, QTableWidgetItem(result_dataset_dict[key]['status']))
                row += 1

   #---------------

    def restart_seqs_homology_search(self, process, result_dataset):
        '''
        Restart a sequences homology search.
        '''

        # initialize the control variable
        OK = True

        # warn that the log window does not have to be closed
        process.write('Do not close this window, please wait!\n')

        # warn that the requirements are being verified
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Checking process requirements ...\n')

        # check the quercusTOA config file
        if OK:
            if not os.path.isfile(genlib.get_app_config_file()):
                process.write(f'*** ERROR: The {genlib.get_app_short_name()} config file does not exist. Please, recreate it.\n')
                OK = False

        # warn that the requirements are OK
        if OK:
            process.write('Process requirements are OK.\n')

        # determine the run directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Determining the run directory ...\n')
            result_dir = self.app_config_dict['Environment parameters']['result_dir']
            current_run_dir = f'{result_dir}/{genlib.get_result_run_subdir()}/{result_dataset}'
            process.write(f'The directory path is {current_run_dir}.\n')

        # set the starter script name
        if OK:
            starter_name = f'{genlib.get_process_search_seqs_homology_code()}-process-starter.sh'

        # submit the starter
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Submitting the starter {starter_name} ...\n')
            command = f'{current_run_dir}/{starter_name} &'
            rc = genlib.run_command(command, process, is_script=True)
            if rc == 0:
                process.write('The script is submitted.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # warn that the log window can be closed
        process.write(f'{genlib.get_separator()}\n')
        process.write('You can close this window now.\n')

        # return the control variable
        return OK

    #---------------

#-------------------------------------------------------------------------------

class FormBrowsetHomologySearch(QWidget):
    '''
    Class used to browse results of a search of sequences homology.
    '''

    #---------------

    def __init__(self, parent):
        '''
        Create a class instance.
        '''

        # save parameters in instance variables
        self.parent = parent

        # call the init method of the parent class
        super().__init__()

        # set the dimensions window
        self.window_height = self.parent.WINDOW_HEIGHT - 100
        self.window_width = self.parent.WINDOW_WIDTH - 50

        # set the head and title
        self.head = 'Browse results of a search of sequences homology'
        self.title = f'{genlib.get_app_short_name()} - {self.head}'

        # get the dictionary of application configuration
        self.app_config_dict = genlib.get_config_dict(genlib.get_app_config_file())

        # build the graphic user interface of the window
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

        # show the window
        self.show()

    #---------------

    def build_gui(self):
        '''
        Build the graphic user interface of the window.
        '''

        # set the width and height of the window
        self.setFixedSize(self.window_width, self.window_height)

        # move the window at center
        rectangle = self.frameGeometry()
        central_point = QGuiApplication.primaryScreen().availableGeometry().center()
        rectangle.moveCenter(central_point)
        self.move(rectangle.topLeft())

        # get font metrics information
        fontmetrics = QFontMetrics(QApplication.font())

        # create and configure "label_head"
        label_head = QLabel(self.head, alignment=Qt.AlignCenter)
        label_head.setStyleSheet('font: bold 14px; color: black; background-color: lightGray; max-height: 30px')

        # create and configure "tablewidget"
        self.tablewidget = QTableWidget()
        self.tablewidget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tablewidget.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.column_name_list = ['Process', 'Result dataset', 'Date', 'Time', 'Status']
        self.tablewidget.setColumnCount(len(self.column_name_list))
        self.tablewidget.setHorizontalHeaderLabels(self.column_name_list)
        self.tablewidget.setColumnWidth(0, 230)
        self.tablewidget.setColumnWidth(1, 280)
        self.tablewidget.setColumnWidth(2, 85)
        self.tablewidget.setColumnWidth(3, 70)
        self.tablewidget.setColumnWidth(4, 90)
        self.tablewidget.verticalHeader().setVisible(True)
        self.tablewidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.tablewidget.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tablewidget.currentCellChanged.connect(self.tablewidget_currentCellChanged)
        self.tablewidget.cellClicked.connect(self.tablewidget_cellClicked)
        self.tablewidget.cellDoubleClicked.connect(self.tablewidget_cellDoubleClicked)

        # create and configure "label_fasta_type"
        label_fasta_type = QLabel()
        label_fasta_type.setText('FASTA type')
        label_fasta_type.setFixedWidth(fontmetrics.width('9'*10))

        # create and configure "lineedit_fasta_type"
        self.lineedit_fasta_type = QLineEdit()
        self.lineedit_fasta_type.editingFinished.connect(self.check_inputs)
        self.lineedit_fasta_type.setDisabled(True)

        # create and configure "label_fasta_file"
        label_fasta_file = QLabel()
        label_fasta_file.setText('FASTA file')
        label_fasta_file.setFixedWidth(fontmetrics.width('9'*20))

        # create and configure "lineedit_fasta_file"
        self.lineedit_fasta_file = QLineEdit()
        self.lineedit_fasta_file.editingFinished.connect(self.check_inputs)
        self.lineedit_fasta_file.setDisabled(True)

        # create and configure "gridlayout_data"
        gridlayout_data = QGridLayout()
        gridlayout_data.setColumnStretch(0, 1)
        gridlayout_data.setColumnStretch(1, 1)
        gridlayout_data.addWidget(self.tablewidget, 0, 0, 1, 11)
        gridlayout_data.addWidget(label_fasta_type, 1, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_fasta_type, 1, 1, 1, 1)
        gridlayout_data.addWidget(label_fasta_file, 2, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_fasta_file, 2, 1, 1, 10)

        # create and configure "groupbox_data"
        groupbox_data = QGroupBox()
        groupbox_data.setObjectName('groupbox_data')
        groupbox_data.setStyleSheet('QGroupBox#groupbox_data {border: 0px;}')
        groupbox_data.setLayout(gridlayout_data)

        # create and configure "pushbutton_refresh"
        self.pushbutton_refresh = QPushButton('Refresh')
        self.pushbutton_refresh.setToolTip('Update the process list.')
        self.pushbutton_refresh.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_refresh.clicked.connect(self.pushbutton_refresh_clicked)

        # create and configure "pushbutton_execute"
        self.pushbutton_execute = QPushButton('Execute')
        self.pushbutton_execute.setToolTip('Restart the process selected.')
        self.pushbutton_execute.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_execute.clicked.connect(self.pushbutton_execute_clicked)

        # create and configure "pushbutton_close"
        pushbutton_close = QPushButton('Close')
        pushbutton_close.setToolTip('Close the window.')
        pushbutton_close.setCursor(QCursor(Qt.PointingHandCursor))
        pushbutton_close.clicked.connect(self.pushbutton_close_clicked)

        # create and configure "gridlayout_buttons"
        gridlayout_buttons = QGridLayout()
        gridlayout_buttons.setColumnStretch(0, 15)
        gridlayout_buttons.setColumnStretch(1, 1)
        gridlayout_buttons.setColumnStretch(2, 1)
        gridlayout_buttons.setColumnStretch(3, 1)
        gridlayout_buttons.addWidget(self.pushbutton_refresh, 0, 1, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(self.pushbutton_execute, 0, 2, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(pushbutton_close, 0, 3, alignment=Qt.AlignCenter)

        # create and configure "groupbox_buttons"
        groupbox_buttons = QGroupBox()
        groupbox_buttons.setObjectName('groupbox_buttons')
        groupbox_buttons.setStyleSheet('QGroupBox#groupbox_buttons {border: 0px;}')
        groupbox_buttons.setLayout(gridlayout_buttons)

        # create and configure "gridlayout_central"
        gridlayout_central = QGridLayout()
        gridlayout_central.setRowStretch(0, 1)
        gridlayout_central.setRowStretch(1, 1)
        gridlayout_central.setRowStretch(2, 10)
        gridlayout_central.setRowStretch(3, 1)
        gridlayout_central.setColumnStretch(0, 1)
        gridlayout_central.addWidget(label_head, 0, 0)
        gridlayout_central.addWidget(QLabel(), 1, 0)
        gridlayout_central.addWidget(groupbox_data, 2, 0)
        gridlayout_central.addWidget(groupbox_buttons, 3, 0)

        # create and configure "groupbox_central"
        groupbox_central = QGroupBox()
        groupbox_central.setLayout(gridlayout_central)

        # create and configure "vboxlayout"
        vboxlayout = QVBoxLayout(self)
        vboxlayout.addWidget(groupbox_central)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load data in "tablewidget"
        self.load_tablewidget()

        # initialize "lineedit_fasta_type"
        self.lineedit_fasta_type.setText('')

        # initialize "lineedit_fasta_file"
        self.lineedit_fasta_file.setText('')

    #---------------

    def check_inputs(self):
        '''
        Check the content of each input and do the actions linked to its value.
        '''

        # initialize the control variable
        OK = True

        # get the list of rows selected
        row_list = []
        for idx in self.tablewidget.selectionModel().selectedIndexes():
            row_list.append(idx.row())
        row_list = list(set(row_list))

        # check "lineedit_fasta_type" when the editing finished
        if not self.lineedit_fasta_type_editing_finished():
            OK = False

        # check "lineedit_fasta_file" when the editing finished
        if not self.lineedit_fasta_file_editing_finished():
            OK = False

        # check all inputs are OK
        if OK:
            self.parent.statusBar().showMessage('')
        else:
            self.parent.statusBar().showMessage('There are one or more inputs without data or with wrong value.')

        # enable "pushbutton_execute"
        if OK and len(row_list) == 1 and self.lineedit_fasta_type.text() != '' and self.lineedit_fasta_file.text() != '':
            self.pushbutton_execute.setEnabled(True)
        else:
            self.pushbutton_execute.setEnabled(False)
            OK = False

        # return the control variable
        return OK

    #---------------

    def tablewidget_currentCellChanged(self, row, _):
        '''
        Perform necessary actions after changing the current "tablewidget" cell.
        '''

        # check if there is a row selected
        if row >= 0:

            # get the dictionary of application configuration
            result_dir = self.app_config_dict['Environment parameters']['result_dir']

            # get the result dataset corresponding to the row clicked
            result_dataset = self.tablewidget.item(row, 1).text()

            # get the parameters file corresponding to the result dataset
            params_file = f'{result_dir}{os.sep}{genlib.get_result_run_subdir()}{os.sep}{result_dataset}{os.sep}{genlib.get_params_file_name()}'
            if sys.platform.startswith('win32'):
                params_file = genlib.wsl_path_2_windows_path(params_file)

            # get the dictionary of parameters
            params_dict = genlib.get_config_dict(params_file)

            # set the FASTA file in "fasta_type"
            fasta_type = params_dict['Search parameters']['fasta_type']
            self.lineedit_fasta_type.setText(fasta_type)

            # set the FASTA file in "fasta_file"
            fasta_file = params_dict['Search parameters']['fasta_file']
            if sys.platform.startswith('win32'):
                fasta_file = genlib.wsl_path_2_windows_path(fasta_file)
            self.lineedit_fasta_file.setText(fasta_file)

        # check the content of inputs
        self.check_inputs()

    #---------------

    def tablewidget_cellClicked(self, _, __):
        '''
        Perform necessary actions after clicking on a "tablewidget" cell.
        '''

        # check the content of inputs
        self.check_inputs()

    #---------------

    def tablewidget_cellDoubleClicked(self, _, __):
        '''
        Perform necessary actions after double clicking on "tablewidget" cell.
        '''

        # check the content of inputs
        OK = self.check_inputs()

        # if inputs are OK, simulate a click on "pushbutton_execute"
        if OK:
            self.pushbutton_execute_clicked()

    #---------------

    def lineedit_fasta_type_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_fasta_type"
        '''

        # initialize the control variable
        OK = True

        # return the control variable
        return OK

    #---------------

    def lineedit_fasta_file_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_fasta_file"
        '''

        # initialize the control variable
        OK = True

        # return the control variable
        return OK

    #---------------

    def pushbutton_refresh_clicked(self):
        '''
        Refresh "tablewidget" and initialize other inputs.
        '''

        # check the content of inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

    #---------------

    def pushbutton_execute_clicked(self):
        '''
        Execute the process.
        '''

        # initialize the control variable
        OK = True

        # get the list of rows selected
        row_list = []
        for idx in self.tablewidget.selectionModel().selectedIndexes():
            row_list.append(idx.row())
        row_list = list(set(row_list))

        # check if there is only a row selected
        if len(row_list) != 1:
            title = f'{genlib.get_app_short_name()} - {self.head}'
            text = 'One row has to be selected.'
            QMessageBox.critical(self, title, text, buttons=QMessageBox.Ok)
            OK = False

        # confirm the process is executed
        if OK:
            text = 'The homology relationships are going to be browsed.\n\nAre you sure to continue?'
            botton = QMessageBox.question(self, self.title, text, buttons=QMessageBox.Yes|QMessageBox.No, defaultButton=QMessageBox.No)
            if botton == QMessageBox.No:
                OK = False

        # execute the process
        if OK:

            # set the process type
            process_type = genlib.get_result_run_subdir()

            # get the result directory
            result_dir = self.app_config_dict['Environment parameters']['result_dir']

            # get the identification of the annotation pipeline dataset
            result_dataset_id = self.tablewidget.item(row_list[0], 1).text()

            # get the file path of the homology relationships
            homology_relationships_file_path = f'{result_dir}{os.sep}{process_type}{os.sep}{result_dataset_id}{os.sep}{genlib.get_homology_relationships_file_name()}'
            if sys.platform.startswith('win32'):
                homology_relationships_file_path = genlib.wsl_path_2_windows_path(homology_relationships_file_path)

            # get the sequences alignment directory corresponding to  the homology relationships
            seq_alignment_dir_path = f'{result_dir}{os.sep}{process_type}{os.sep}{result_dataset_id}{os.sep}seqs-alignments'
            if sys.platform.startswith('win32'):
                seq_alignment_dir_path = genlib.wsl_path_2_windows_path(seq_alignment_dir_path)

            # get homology relationships
            QApplication.setOverrideCursor(Qt.WaitCursor)
            (homology_relationships_dict, data_list, data_dict, selection_data, window_height, window_width, explanatory_text) = self.get_homology_relationships(homology_relationships_file_path)
            QGuiApplication.restoreOverrideCursor()

            # show homology relationships
            head = f'Homology relationships file {homology_relationships_file_path}'
            data_table = dialogs.DialogHomologyRelationships(self, head, window_height, window_width, data_list, data_dict, selection_data, homology_relationships_dict, sorted(homology_relationships_dict.keys()), seq_alignment_dir_path)
            data_table.exec()

        # close the windows
        # -- if OK:
        # --     self.pushbutton_close_clicked()

    #---------------

    def pushbutton_close_clicked(self):
        '''
        Close the window.
        '''

        self.parent.current_subwindow = None
        self.close()
        self.parent.set_background_image()

    #---------------

    def load_tablewidget(self):
        '''
        Load data in "tablewidget".
        '''

        # get the result directory
        result_dir = self.app_config_dict['Environment parameters']['result_dir']

        # set the type, name and code of the annotation pipeline datasets
        process_type = genlib.get_result_run_subdir()
        process_name = genlib.get_process_search_seqs_homology_name()
        process_code = genlib.get_process_id(process_name)

        # get the process dictionary
        process_dict = genlib.get_process_dict()

        # get the log directory
        log_dir = f'{result_dir}/{process_type}'
        if sys.platform.startswith('win32'):
            log_dir = genlib.wsl_path_2_windows_path(log_dir)

        # set the command to get the result datasets of annotation pipeline in the log directory
        command = ''
        if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
            if process_name == 'all':
                command = f'ls -d {log_dir}/*  | xargs -n 1 basename'
            else:
                command = f'ls -d {log_dir}/{process_code}-*  | xargs -n 1 basename'
        elif sys.platform.startswith('win32'):
            log_dir = log_dir.replace('/', '\\')
            if process_name == 'all':
                command = f'dir /a:d /b {log_dir}\\*'
            else:
                command = f'dir /a:d /b {log_dir}\\{process_code}-*'

        # run the command to get the result datasets of enrichment analysis in the log directory
        output = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=False)

        # initialize the result dataset dictionary
        result_dataset_dict = {}

        # build the result dataset dictionary
        for line in output.stdout.split('\n'):
            if line != '':

                # get data
                result_dataset_id = line.strip()
                try:
                    pattern = r'^(.+)\-(.+)\-(.+)$'
                    mo = re.search(pattern, result_dataset_id)
                    process_code = mo.group(1).strip()
                    process_name = process_dict[process_code]['name']
                    yymmdd = mo.group(2)
                    hhmmss = mo.group(3)
                    date = f'20{yymmdd[:2]}-{yymmdd[2:4]}-{yymmdd[4:]}'
                    time = f'{hhmmss[:2]}:{hhmmss[2:4]}:{hhmmss[4:]}'
                except:    # pylint: disable=bare-except
                    process_name = 'unknown process'
                    date = '0000-00-00'
                    time = '00:00:00'

                # determine the status
                status_ok = os.path.isfile(genlib.get_status_ok(os.path.join(log_dir, result_dataset_id)))
                status_wrong = os.path.isfile(genlib.get_status_wrong(os.path.join(log_dir, result_dataset_id)))
                status = ''
                if status_ok and not status_wrong:
                    status = 'OK'
                elif not status_ok and status_wrong:
                    status = 'wrong'
                elif not status_ok and not status_wrong:
                    status = 'not finished'
                elif status_ok and status_wrong:
                    status = 'undetermined'

                # insert data in the dictionary when the the status is OK
                if status == 'OK':
                    key = f'{process_name}-{result_dataset_id}'
                    result_dataset_dict[key] = {'process': process_name, 'result_dataset_id': result_dataset_id, 'date': date, 'time': time, 'status': status}

        # initialize "tablewidget"
        self.tablewidget.clearContents()

        # set the rows number of "tablewidget"
        self.tablewidget.setRowCount(len(result_dataset_dict))

        # load data in "tablewidget" for the OK result datasets of annotation pipeline
        if not result_dataset_dict:
            text = 'There is no run ended OK.'
            QMessageBox.warning(self, self.title, text, buttons=QMessageBox.Ok)
        else:
            row = 0
            for key in sorted(result_dataset_dict.keys()):
                self.tablewidget.setItem(row, 0, QTableWidgetItem(result_dataset_dict[key]['process']))
                self.tablewidget.setItem(row, 1, QTableWidgetItem(result_dataset_dict[key]['result_dataset_id']))
                self.tablewidget.setItem(row, 2, QTableWidgetItem(result_dataset_dict[key]['date']))
                self.tablewidget.setItem(row, 3, QTableWidgetItem(result_dataset_dict[key]['time']))
                self.tablewidget.setItem(row, 4, QTableWidgetItem(result_dataset_dict[key]['status']))
                row += 1

   #---------------

    @staticmethod
    def get_homology_relationships(homology_relationships_file):
        '''
        Get homology relationships.
        '''

        # initialize the homology relationships dictionary
        homology_relatioships_dict = {}

        # open the homology relationships file
        if homology_relationships_file.endswith('.gz'):
            try:
                homology_relationships_file_id = gzip.open(homology_relationships_file, mode='rt', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                raise genlib.ProgramException(e, 'F002', homology_relationships_file)
        else:
            try:
                homology_relationships_file_id = open(homology_relationships_file, mode='r', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                raise genlib.ProgramException(e, 'F001', homology_relationships_file)

        # initialize the record counter
        record_counter = 0

        # initialize the header record control
        header_record = True

        # read the first record
        record = homology_relationships_file_id.readline()

        # while there are records
        while record != '':

            # add 1 to the record counter
            record_counter += 1

            # process the header record
            if header_record:
                header_record = False

            # process data records
            else:

                # extract data
                # record format: Sequence id;Species name;Homologous gene id;Homologous protein isoforms
                data_list = []
                begin = 0
                for end in [i for i, chr in enumerate(record) if chr == ';']:
                    data_list.append(record[begin:end])
                    begin = end + 1
                data_list.append(record[begin:].strip('\n'))
                try:
                    seq_id = data_list[0]
                    species_name = genlib.get_quercus_species_name(data_list[1])
                    homologous_gene_id = data_list[2]
                    homologous_protein_isoforms = data_list[3]
                except Exception as e:
                    raise genlib.ProgramException(e, 'F006', os.path.basename(homology_relationships_file), record_counter)

                # add data to the dictionary
                key = f'{seq_id}-{species_name}'
                homology_relatioships_dict[key] = {'seq_id': seq_id, 'species_name': species_name, 'homologous_gene_id': homologous_gene_id, 'homologous_protein_isoforms': homologous_protein_isoforms}

            # read the next record
            record = homology_relationships_file_id.readline()

        # build the data list
        data_list = ['seq_id', 'species_name', 'homologous_gene_id', 'homologous_protein_isoforms']

        # build the data dictionary
        data_dict = {}
        data_dict['seq_id'] = {'text': 'Sequence id', 'width': 150, 'alignment': 'left'}
        data_dict['species_name'] = {'text': 'Species name', 'width': 150, 'alignment': 'left'}
        data_dict['homologous_gene_id'] = {'text': 'Homologous gene id', 'width': 180, 'alignment': 'left'}
        data_dict['homologous_protein_isoforms'] = {'text': 'Homologous protein isoforms', 'width': 600, 'alignment': 'left'}

        # set the selection data
        selection_data = 'seq_id'

        # set the window height and width
        window_height = 600
        window_width = 1200

        # set the explanatory text
        explanatory_text = ''

        # return data
        return homology_relatioships_dict, data_list, data_dict, selection_data, window_height, window_width, explanatory_text

    #---------------

#-------------------------------------------------------------------------------

class FormGetProteinIdsHomology(QWidget):
    '''
    Class used to get homology data of protein identifications (accession numbers).
    '''

    #---------------

    def __init__(self, parent):
        '''
        Create a class instance.
        '''

        # save parameters in instance variables
        self.parent = parent

        # call the init method of the parent class
        super().__init__()

        # set the dimensions window
        self.window_height = self.parent.WINDOW_HEIGHT - 100
        self.window_width = self.parent.WINDOW_WIDTH - 50

        # set the head and title
        self.head = 'Get homology data of protein identifications (accession numbers)'
        self.title = f'{genlib.get_app_short_name()} - {self.head}'

        # get the dictionary of application configuration
        self.app_config_dict = genlib.get_config_dict(genlib.get_app_config_file())

        # connect to the sequences database
        sequences_database = self.app_config_dict['quercusTOA-app database']['sequences_db_path']
        if sys.platform.startswith('win32'):
            sequences_database = genlib.wsl_path_2_windows_path(sequences_database)
        self.sequences_database_conn = sqllib.connect_database(sequences_database)

        # connect to the comparative genomics database
        comparative_genomics_database = self.app_config_dict['quercusTOA-app database']['comparative_genomics_db_path']
        if sys.platform.startswith('win32'):
            comparative_genomics_database = genlib.wsl_path_2_windows_path(comparative_genomics_database)
        self.comparative_genomics_database_conn = sqllib.connect_database(comparative_genomics_database)

        # build the graphic user interface of the window
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

        # show the window
        self.show()

    #---------------

    def build_gui(self):
        '''
        Build the graphic user interface of the window.
        '''

        # set the width and height of the window
        self.setFixedSize(self.window_width, self.window_height)

        # move the window at center
        rectangle = self.frameGeometry()
        central_point = QGuiApplication.primaryScreen().availableGeometry().center()
        rectangle.moveCenter(central_point)
        self.move(rectangle.topLeft())

        # get font metrics information
        fontmetrics = QFontMetrics(QApplication.font())

        # create and configure "label_head"
        label_head = QLabel(self.head, alignment=Qt.AlignCenter)
        label_head.setStyleSheet('font: bold 14px; color: black; background-color: lightGray; max-height: 30px')

        # create and configure "label_protein_ids"
        label_protein_ids = QLabel()
        label_protein_ids.setText('Protein identifications (accession numbers)\nseparated by commas')
        label_protein_ids.setFixedWidth(fontmetrics.width('9'*35))

        # create and configure "plaintextedit_protein_ids_list"
        self.plaintextedit_protein_ids_list = QPlainTextEdit()
        self.plaintextedit_protein_ids_list.setFixedWidth(fontmetrics.width('9'*75))
        self.plaintextedit_protein_ids_list.setFixedHeight(fontmetrics.lineSpacing() * 3 + 12)
        self.plaintextedit_protein_ids_list.textChanged.connect(self.check_inputs)

        # create and configure "gridlayout_data"
        gridlayout_data = QGridLayout()
        gridlayout_data.setColumnStretch(0,1)
        gridlayout_data.setColumnStretch(1,1)
        gridlayout_data.addWidget(label_protein_ids, 0, 0, 1, 1)
        gridlayout_data.addWidget(self.plaintextedit_protein_ids_list, 0, 1, 1, 1, alignment=Qt.AlignLeft)

        # create and configure "groupbox_data"
        groupbox_data = QGroupBox()
        groupbox_data.setObjectName('groupbox_data')
        groupbox_data.setStyleSheet('QGroupBox#groupbox_data {border: 0px;}')
        groupbox_data.setLayout(gridlayout_data)

        # create and configure "pushbutton_execute"
        self.pushbutton_execute = QPushButton('Execute')
        self.pushbutton_execute.setToolTip('Execute the search of homology data.')
        self.pushbutton_execute.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_execute.clicked.connect(self.pushbutton_execute_clicked)

        # create and configure "pushbutton_close"
        pushbutton_close = QPushButton('Close')
        pushbutton_close.setToolTip('Cancel the the search of homology data and close the window.')
        pushbutton_close.setCursor(QCursor(Qt.PointingHandCursor))
        pushbutton_close.clicked.connect(self.pushbutton_close_clicked)

        # create and configure "gridlayout_buttons"
        gridlayout_buttons = QGridLayout()
        gridlayout_buttons.setColumnStretch(0, 10)
        gridlayout_buttons.setColumnStretch(1, 1)
        gridlayout_buttons.setColumnStretch(2, 1)
        gridlayout_buttons.addWidget(self.pushbutton_execute, 0, 1, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(pushbutton_close, 0, 2, alignment=Qt.AlignCenter)

        # create and configure "groupbox_buttons"
        groupbox_buttons = QGroupBox()
        groupbox_buttons.setObjectName('groupbox_buttons')
        groupbox_buttons.setStyleSheet('QGroupBox#groupbox_buttons {border: 0px;}')
        groupbox_buttons.setLayout(gridlayout_buttons)

        # create and configure "gridlayout_central"
        gridlayout_central = QGridLayout()
        gridlayout_central.setRowStretch(0, 1)
        gridlayout_central.setRowStretch(1, 1)
        gridlayout_central.setRowStretch(2, 1)
        gridlayout_central.setRowStretch(3, 1)
        gridlayout_central.setRowStretch(4, 1)
        gridlayout_central.setColumnStretch(0, 0)
        gridlayout_central.setColumnStretch(1, 1)
        gridlayout_central.setColumnStretch(2, 0)
        gridlayout_central.addWidget(label_head, 0, 1)
        gridlayout_central.addWidget(QLabel(), 1, 1)
        gridlayout_central.addWidget(groupbox_data, 2, 1)
        gridlayout_central.addWidget(QLabel(), 3, 1)
        gridlayout_central.addWidget(groupbox_buttons, 4, 1)

        # create and configure "groupbox_central"
        groupbox_central = QGroupBox()
        groupbox_central.setLayout(gridlayout_central)

        # create and configure "vboxlayout"
        vboxlayout = QVBoxLayout(self)
        vboxlayout.addWidget(groupbox_central)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # set initial value in "plaintextedit_protein_ids_list"
        self.plaintextedit_protein_ids_list.setPlaceholderText('')

    #---------------

    def check_inputs(self):
        '''
        Check the content of each input and do the actions linked to its value.
        '''

        # initialize the control variable
        OK = True

        # check "plaintextedit_protein_ids_list" when the editing finished
        if not self.plaintextedit_protein_ids_list_editing_finished():
            OK = False

        # enable "pushbutton_execute"
        if OK and self.plaintextedit_protein_ids_list.toPlainText() != '':
            self.pushbutton_execute.setEnabled(True)
        else:
            self.pushbutton_execute.setEnabled(False)

        # return the control variable
        return OK

    #---------------

    def plaintextedit_protein_ids_list_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "plaintextedit_protein_ids_list"
        '''

        # initialize the control variable
        OK = True

        # chek if "plaintextedit_protein_ids_list" is empty
        if self.plaintextedit_protein_ids_list.toPlainText() == '':
            OK = False
            self.plaintextedit_protein_ids_list.setStyleSheet('background-color: white')

        # return the control variable
        return OK

    #---------------

    def pushbutton_execute_clicked(self):
        '''
        Execute the process.
        '''

        # initialize the control variable
        OK = True

        # confirm the process is executed
        text = ''
        if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
            text = 'The process to search the homology data is going to be run.\n\nIt may take some time due to the online alignments.\n\nAre you sure to continue?'
        elif sys.platform.startswith('win32'):
            text = 'The process to search the homology data is going to be run.\n\nAre you sure to continue?'
        botton = QMessageBox.question(self, self.title, text, buttons=QMessageBox.Yes|QMessageBox.No, defaultButton=QMessageBox.No)
        if botton == QMessageBox.No:
            OK = False

        # execute the process
        if OK:

            # get the protein_id
            reference_protein_ids_list = [x.strip() for x in self.plaintextedit_protein_ids_list.toPlainText().split(',')]

            # get homology data
            QApplication.setOverrideCursor(Qt.WaitCursor)
            (homology_relationships_dict, data_list, data_dict, selection_data, window_height, window_width, explanatory_text) = self.get_homology_relationships(reference_protein_ids_list)
            QGuiApplication.restoreOverrideCursor()

            # check if there are item in the homology data dictionary
            if homology_relationships_dict:

                if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):

                    # build alignments of the protein and gene FASTA sequences corresponding to homology relationships
                    QApplication.setOverrideCursor(Qt.WaitCursor)
                    OK = self.build_alignments(homology_relationships_dict)
                    QGuiApplication.restoreOverrideCursor()

                    # get the application directory
                    app_dir = self.app_config_dict['Environment parameters']['app_dir']

                    # set the directory where the FASTA sequence and alignment files are located
                    seq_alignment_dir_path = f'{app_dir}{os.sep}{genlib.get_temp_dir()}'

                    # show dialog with homology relationships
                    head = 'Homology relationships'
                    data_table = dialogs.DialogHomologyRelationships(self, head, window_height, window_width, data_list, data_dict, selection_data, homology_relationships_dict, sorted(homology_relationships_dict.keys()), seq_alignment_dir_path)
                    data_table.exec()

                elif sys.platform.startswith('win32'):

                    # show dialog with homology relationships
                    head = 'Homology relationships'
                    data_table = dialogs.DialogDataTableWithSelections(self, head, window_height, window_width, data_list, data_dict, selection_data, homology_relationships_dict, sorted(homology_relationships_dict.keys()), explanatory_text, '')
                    data_table.exec()

            else:

                # show message
                title = f'{genlib.get_app_short_name()} - Warning'
                text = f'No homology relationships have been found for these protein identifications (accession numbers) in {genlib.get_db_name()}.'
                QMessageBox.warning(self, title, text, buttons=QMessageBox.Ok)

        # close the windows
        # -- if OK:
        # --     self.pushbutton_close_clicked()

    #---------------

    def pushbutton_close_clicked(self):
        '''
        Close the window.
        '''

        self.parent.current_subwindow = None
        self.close()
        self.parent.set_background_image()

    #---------------

    def get_homology_relationships(self, reference_protein_ids_list):
        '''
        Get homology relationships corresponding to the list of reference protein identifications.
        '''

        # initialize the homology relationships dictionary
        homology_relationships_dict = genlib.NestedDefaultDict()

        # for each identification in the list of reference protein identifications
        for reference_protein_id in sorted(reference_protein_ids_list):

            # get the list of protein isoforms corresponding to the reference protein identification
            mmseqs2_protein_isoforms_list = sqllib.get_mmseqs2_protein_isoforms_list(self.comparative_genomics_database_conn, '', [reference_protein_id])

            # build the list of protein isoforms identification corresponding to the reference protein identification
            reference_protein_isoform_ids_list = []
            for _, protein_isoform_data in enumerate(mmseqs2_protein_isoforms_list):
                reference_protein_isoform_ids_list.append(protein_isoform_data['protein_id'])

            # check if there are items in the list of protein isoforms corresponding to the reference protein identification
            if reference_protein_isoform_ids_list:

                # add data of reference protein to the homology relationships dictionary
                reference_protein_isoform_ids_list_reviewed =  [x for x in reference_protein_isoform_ids_list if x != reference_protein_id]
                if reference_protein_isoform_ids_list_reviewed:
                    species_id = mmseqs2_protein_isoforms_list[0]['species_id']
                    species_name = genlib.get_quercus_species_name(species_id)
                    gene_id = mmseqs2_protein_isoforms_list[0]['gene_id']
                    protein_isoform_ids = '|'.join(sorted(reference_protein_isoform_ids_list_reviewed))
                    homology_relationships_dict[f'{reference_protein_id}-{species_id}-{gene_id}'] = {'reference_protein_id': reference_protein_id, 'species_id': species_id, 'species_name': species_name, 'gene_id':gene_id, 'protein_isoform_ids': protein_isoform_ids}

                # for each identification in the list of protein isoforms corresponding to the reference protein identification
                for reference_protein_isoform_id in reference_protein_isoform_ids_list:

                    # get data of the orthologous proteins
                    orthologous_protein_data_list = sqllib.get_orthologous_protein_data_list(self.comparative_genomics_database_conn, reference_protein_isoform_id)

                    # check if there are data of orthologous proteins:
                    if orthologous_protein_data_list:

                        # cluster the list of target orthologous proteins by species identification and gene identification
                        clustered_homologous_proteins_dict = defaultdict(list)
                        for item in orthologous_protein_data_list:
                            key = (item['target_species_id'], item['gene_id'])
                            clustered_homologous_proteins_dict[key].append(item['target_protein_id'])

                        # build the orthologous data list
                        orthologous_data_list = [
                            {
                                'species_id': species_id,
                                'gene_id': gene_id,
                                'protein_isoform_ids': '|'.join(sorted(protein_ids_list))
                            }
                            for (species_id, gene_id), protein_ids_list in clustered_homologous_proteins_dict.items()
                        ]

                        # add orthologous data to the homology relationships dictionary
                        for _, data in enumerate(orthologous_data_list):
                            species_id = data['species_id']
                            species_name = genlib.get_quercus_species_name(species_id)
                            gene_id = data['gene_id']
                            protein_isoform_ids = data['protein_isoform_ids']
                            homology_relationships_dict[f'{reference_protein_id}-{species_id}-{gene_id}'] = {'reference_protein_id': reference_protein_id, 'species_id': species_id, 'species_name': species_name, 'gene_id':gene_id, 'protein_isoform_ids': protein_isoform_ids}

        # build the data list
        data_list = ['reference_protein_id', 'species_name', 'gene_id', 'protein_isoform_ids']

        # build the data dictionary
        data_dict = {}
        data_dict['reference_protein_id'] = {'text': 'Protein id', 'width': 150, 'alignment': 'left'}
        data_dict['species_name'] = {'text': 'Species name', 'width': 150, 'alignment': 'left'}
        data_dict['gene_id'] = {'text': 'Homologous gene id', 'width': 180, 'alignment': 'left'}
        data_dict['protein_isoform_ids'] = {'text': 'Homologous protein isoforms', 'width': 600, 'alignment': 'left'}

        # set the selection data
        selection_data = 'reference_protein_id'

        # set the window height and width
        window_height = 600
        window_width = 1200

        # set the explanatory text
        explanatory_text = ''

        # return data
        return homology_relationships_dict, data_list, data_dict, selection_data, window_height, window_width, explanatory_text

    #---------------

    def build_alignments(self, homology_relationships_dict):
        '''
        Build alignments of the protein and gene FASTA sequences corresponding to homology relationships.
        '''

        # initialize the control variable
        OK = True

        # set the homology relationships file
        homology_relationships_file = f'{genlib.get_temp_dir()}{os.sep}{genlib.get_homology_relationships_file_name()}'

        # save data in the homology relationships file
        with open(homology_relationships_file, mode='w', encoding='iso-8859-1') as homology_relationships_file_id:
            for key in sorted(homology_relationships_dict):
                reference_protein_id = homology_relationships_dict[key]['reference_protein_id']
                species_id = homology_relationships_dict[key]['species_id']
                gene_id = homology_relationships_dict[key]['gene_id']
                protein_isoform_ids = homology_relationships_dict[key]['protein_isoform_ids']
                homology_relationships_file_id.write(f'{reference_protein_id};{species_id};{gene_id};{protein_isoform_ids}\n')

        # open the homology relationships file
        try:
            homology_relationships_file_id = open(homology_relationships_file, mode='r', encoding='iso-8859-1')
        except Exception as e:
            raise genlib.ProgramException(e, 'F001', homology_relationships_file)

        # initialize the counter of records corresponding to the homology relationships file
        homology_relationships_record_counter = 0

        # read the first record of the homology relationships file
        (homology_relationships_record, _, homology_relationships_data_dict) = genlib.read_homology_relationships_record(homology_relationships_file, homology_relationships_file_id, homology_relationships_record_counter)

        # while there are records in the homology relationships file
        while homology_relationships_record != '' and OK:

            # save the sequence identification
            w_seq_id = homology_relationships_data_dict['seq_id']

            # set the file path of homologous protein FASTA sequences corresponding to the current sequence identification
            protein_seq_fasta_file = f'{genlib.get_temp_dir()}{os.sep}{w_seq_id}-homologous-proteins.fasta'

            # open the protein sequence FASTA file
            try:
                protein_seq_fasta_file_id = open(protein_seq_fasta_file, mode='w', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                raise genlib.ProgramException(e, 'F003', protein_seq_fasta_file)

            # set the file path of homologous gene FASTA sequences corresponding to the current sequence identification
            gene_seq_fasta_file = f'{genlib.get_temp_dir()}{os.sep}{w_seq_id}-homologous-genes.fasta'

            # open the gene sequence FASTA file
            try:
                gene_seq_fasta_file_id = open(gene_seq_fasta_file, mode='w', encoding='iso-8859-1', newline='\n')
            except Exception as e:
                raise genlib.ProgramException(e, 'F003', gene_seq_fasta_file)

            # get the sequence of the protein identification
            protein_id_seq_dict = sqllib.get_protein_seq_dict(self.sequences_database_conn, w_seq_id)

            # write the sequence of the protein identification in the the protein sequence FASTA file
            protein_seq_fasta_file_id.write(f'>{w_seq_id}[{protein_id_seq_dict['species_id']}]\n')
            protein_seq_fasta_file_id.write(f'{protein_id_seq_dict['seq']}\n')

            # get the gene identification and species corresponding to que protein identification
            protein_data_dict = sqllib.get_mmseqs2_protein_data_dict(self.comparative_genomics_database_conn, w_seq_id)
            protein_gene_id = protein_data_dict['gene_id']
            species_id = protein_data_dict['species_id']

            # get the gene sequence data of protein identification
            gene_seq_dict = sqllib.get_gene_seq_dict(self.sequences_database_conn, protein_gene_id)
            gene_seq = gene_seq_dict['seq']

            # write the gene sequence of the protein identification in the  the protein FASTA sequence file
            gene_seq_fasta_file_id.write(f'>{protein_gene_id}[{species_id}]\n')
            gene_seq_fasta_file_id.write(f'{gene_seq}\n')

            # while there are records in the homology relationships file and same sequence identification
            while homology_relationships_record != '' and homology_relationships_data_dict['seq_id'] ==  w_seq_id:

                # add 1 to the counter of records corresponding to the homology relationships file
                homology_relationships_record_counter += 1

                # set the protein isoform identification list
                protein_isoform_id_list = homology_relationships_data_dict['protein_isoform_ids'].split('|')

                # for each protein isoform identification
                for protein_isoform_id in sorted(protein_isoform_id_list):

                    # get the protein isoform sequence data
                    protein_isoform_seq_dict = sqllib.get_protein_seq_dict(self.sequences_database_conn, protein_isoform_id)
                    protein_isoform_seq = protein_isoform_seq_dict['seq']

                    # write the protein isoform sequence in the protein sequence FASTA file
                    protein_seq_fasta_file_id.write(f'>{protein_isoform_id}[{homology_relationships_data_dict['species_id']}]\n')
                    protein_seq_fasta_file_id.write(f'{protein_isoform_seq}\n')

                    # print(f"homology_relationships_data_dict['gene_id']: {homology_relationships_data_dict['gene_id']}")

                    # get the gene sequence data of protein isoform
                    gene_seq_dict = sqllib.get_gene_seq_dict(self.sequences_database_conn, homology_relationships_data_dict['gene_id'])
                    gene_isoform_seq = gene_seq_dict['seq']

                    # print(f"gene_isoform_seq: {gene_isoform_seq}")

                    # write the gene sequence  of the protein isoform in the  the protein FASTA sequence file
                    gene_seq_fasta_file_id.write(f'>{homology_relationships_data_dict['gene_id']}[{homology_relationships_data_dict['species_id']}]\n')
                    gene_seq_fasta_file_id.write(f'{gene_isoform_seq}\n')

                # print the record counters
                genlib.Message.print('verbose', f'\rHolology relationships records: {homology_relationships_record_counter:5d}')

                # read the next record of the homology relationships file
                (homology_relationships_record, _, homology_relationships_data_dict) = genlib.read_homology_relationships_record(homology_relationships_file, homology_relationships_file_id, homology_relationships_record_counter)

            # close sequence FASTA files of the current sequence identification
            protein_seq_fasta_file_id.close()
            gene_seq_fasta_file_id.close()

            # get the sequence number in the protein FASTA sequence file
            with open(protein_seq_fasta_file, 'r', encoding='iso-8859-1') as protein_seq_fasta_file_id:
                record_count = sum(1 for _ in protein_seq_fasta_file_id)
            protein_seq_number = record_count / 2

            # set the protein alignment file
            protein_alignment_file = f'{genlib.get_temp_dir()}{os.sep}{w_seq_id}-homologous-proteins.fasta.aln'

            # run alignment of their homologous protein isoforms
            if protein_seq_number > 1:
                with open(protein_alignment_file, mode='w', encoding='iso-8859-1') as alignment_file_id:
                    result = subprocess.run(['conda', 'run', '-n', genlib.get_mafft_environment(),
                                            'mafft', '--auto', '--anysymbol', protein_seq_fasta_file],
                                            stdout=alignment_file_id, stderr=subprocess.PIPE, text=True, check=True)
                if result.returncode != 0:
                    OK = False
                    title = f'{genlib.get_app_short_name()} - MAFFT alignment'
                    text = f' error:\n\n{result.stderr}.'
                    QMessageBox.critical(self, title, text, buttons=QMessageBox.Ok)
            else:
                shutil.copy(protein_seq_fasta_file, protein_alignment_file)

            # plot the alignment of homologous protein isoform using pyMSAviz
            if OK:
                protein_alignment_plot_file = f'{genlib.get_temp_dir()}{os.sep}{w_seq_id}-homologous-proteins.fasta.aln.pdf'
                protein_alignment_plot = MsaViz(protein_alignment_file, format='fasta', wrap_length=80, sort=False, color_scheme='Flower', show_count=True, show_consensus=False)
                protein_alignment_plot.savefig(protein_alignment_plot_file)

            # generate the guide tree represents the clustering of homologous protein isoform sequences in Newick format when there is more than one FASTA sequence
            if OK and protein_seq_number > 1:
                result = subprocess.run(['conda', 'run', '-n', genlib.get_mafft_environment(),
                                        'mafft', '--treeout', protein_seq_fasta_file],
                                        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, check=True)
                if result.returncode != 0:
                    OK = False
                    title = f'{genlib.get_app_short_name()} - MAFFT guide tree generation'
                    text = f' error:\n\n{result.stderr}.'
                    QMessageBox.critical(self, title, text, buttons=QMessageBox.Ok)

            # plot the guide tree when there is more than one FASTA sequence
            if OK and protein_seq_number > 1:
                tree = Phylo.read(f'{protein_seq_fasta_file}.tree', 'newick')
                Phylo.draw(tree, do_show=False)
                tree_plot_file = f'{protein_seq_fasta_file}.tree.pdf'
                plt.savefig(tree_plot_file)
                plt.close()

            # get the sequence number in the homologous protein isoform genes
            if OK:
                with open(gene_seq_fasta_file, 'r', encoding='iso-8859-1') as gene_seq_fasta_file_id:
                    record_count = sum(1 for _ in gene_seq_fasta_file_id)
                gene_seq_number = record_count / 2

            # run alignment of their homologous protein isoform genes
            if OK:
                if gene_seq_number > 1:
                    gene_alignment_file = f'{genlib.get_temp_dir()}{os.sep}{w_seq_id}-homologous-genes.fasta.aln'
                    with open(gene_alignment_file, mode='w', encoding='iso-8859-1') as alignment_file_id:
                        result = subprocess.run(['conda', 'run', '-n', genlib.get_mafft_environment(),
                                                    'mafft', '--auto', '--anysymbol', gene_seq_fasta_file],
                                                    stdout=alignment_file_id, stderr=subprocess.PIPE, text=True, check=True)
                    if result.returncode != 0:
                        OK = False
                        title = f'{genlib.get_app_short_name()} - MAFFT alignment'
                        text = f' error:\n\n{result.stderr}.'
                        QMessageBox.critical(self, title, text, buttons=QMessageBox.Ok)
                else:
                    shutil.copy(gene_seq_fasta_file, gene_alignment_file)

            # plot the alignment  using pyMSAviz when there is more than one FASTA sequence
            if OK and gene_seq_number > 1:
                gene_alignment_plot_file = f'{genlib.get_temp_dir()}{os.sep}{w_seq_id}-homologous-genes.fasta.aln.pdf'
                gene_alignment_plot = MsaViz(gene_alignment_file, format='fasta', wrap_length=80, sort=False, color_scheme='Flower', show_count=True, show_consensus=False)
                gene_alignment_plot.savefig(gene_alignment_plot_file)

        # close the homology relationships file
        homology_relationships_file_id.close()

        # return the control variable
        return OK

    #---------------

#-------------------------------------------------------------------------------

class FormProcessGffFileBuilding(QWidget):
    '''
    Class used to create a GFF file of a target species genome using a reference species data.
    '''

    #---------------

    def __init__(self, parent):
        '''
        Create a class instance.
        '''

        # save parameters in instance variables
        self.parent = parent

        # call the init method of the parent class
        super().__init__()

        # set the dimensions window
        self.window_height = self.parent.WINDOW_HEIGHT - 100
        self.window_width = self.parent.WINDOW_WIDTH - 50

        # set the head and title
        self.head = 'Create a GFF file of a target species genome using a reference species data'
        self.title = f'{genlib.get_app_short_name()} - {self.head}'

        # get the dictionary of application configuration
        self.app_config_dict = genlib.get_config_dict(genlib.get_app_config_file())

        # build the graphic user interface of the window
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

        # show the window
        self.show()

    #---------------

    def build_gui(self):
        '''
        Build the graphic user interface of the window.
        '''

        # set the width and height of the window
        self.setFixedSize(self.window_width, self.window_height)

        # move the window at center
        rectangle = self.frameGeometry()
        central_point = QGuiApplication.primaryScreen().availableGeometry().center()
        rectangle.moveCenter(central_point)
        self.move(rectangle.topLeft())

        # get font metrics information
        fontmetrics = QFontMetrics(QApplication.font())

        # create and configure "label_head"
        label_head = QLabel(self.head, alignment=Qt.AlignCenter)
        label_head.setStyleSheet('font: bold 14px; color: black; background-color: lightGray; max-height: 30px')

        # create and configure "label_threads"
        label_threads = QLabel()
        label_threads.setText('Threads')
        label_threads.setFixedWidth(fontmetrics.width('9'*10))

        # create and configure "lineedit_threads"
        self.lineedit_threads  = QLineEdit()
        self.lineedit_threads.setFixedWidth(fontmetrics.width('9'*6))
        self.lineedit_threads.editingFinished.connect(self.check_inputs)

        # create and configure "label_reference_species_name"
        label_reference_species_name = QLabel()
        label_reference_species_name.setText('Reference species')
        label_reference_species_name.setFixedWidth(fontmetrics.width('9'*15))

        # create and configure "combobox_reference_species_name"
        self.combobox_reference_species_name = QComboBox()
        self.combobox_reference_species_name.currentIndexChanged.connect(self.check_inputs)
        self.combobox_reference_species_name.setFixedWidth(fontmetrics.width('9'*22))

        # create and configure "label_target_genome_path"
        label_target_genome_path = QLabel()
        label_target_genome_path.setText('Genome FASTA file')

        # create and configure "lineedit_target_genome_path"
        self.lineedit_target_genome_path = QLineEdit()
        self.lineedit_target_genome_path.editingFinished.connect(self.check_inputs)

        # create and configure "pushbutton_search_target_genome_path"
        pushbutton_search_target_genome_path = QPushButton('Search ...')
        pushbutton_search_target_genome_path.setToolTip('Search and select the genome FASTA file.')
        pushbutton_search_target_genome_path.setCursor(QCursor(Qt.PointingHandCursor))
        pushbutton_search_target_genome_path.clicked.connect(self.pushbutton_search_target_genome_path_clicked)

        # create and configure "empty"
        label_empty = QLabel()
        label_empty.setFixedWidth(fontmetrics.width('9'*3))

        # create and configure "gridlayout_data"
        gridlayout_data = QGridLayout()
        gridlayout_data.setRowMinimumHeight(0, 40)
        gridlayout_data.setRowMinimumHeight(1, 40)
        gridlayout_data.setColumnStretch(0,1)
        gridlayout_data.setColumnStretch(1,1)
        gridlayout_data.setColumnStretch(2,1)
        gridlayout_data.setColumnStretch(3,1)
        gridlayout_data.setColumnStretch(4,15)
        gridlayout_data.setColumnStretch(5,1)
        gridlayout_data.addWidget(label_threads, 0, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_threads, 0, 1, 1, 1, alignment=Qt.AlignLeft)
        gridlayout_data.addWidget(label_empty, 0, 2, 1, 1)
        gridlayout_data.addWidget(label_reference_species_name, 0, 3, 1, 1)
        gridlayout_data.addWidget(self.combobox_reference_species_name, 0, 4, 1, 1, alignment=Qt.AlignLeft)
        gridlayout_data.addWidget(label_target_genome_path, 1, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_target_genome_path, 1, 1, 1, 4)
        gridlayout_data.addWidget(pushbutton_search_target_genome_path, 1, 5, 1, 1)

        # create and configure "groupbox_data"
        groupbox_data = QGroupBox()
        groupbox_data.setObjectName('groupbox_data')
        groupbox_data.setStyleSheet('QGroupBox#groupbox_data {border: 0px;}')
        groupbox_data.setLayout(gridlayout_data)

        # create and configure "pushbutton_execute"
        self.pushbutton_execute = QPushButton('Execute')
        self.pushbutton_execute.setToolTip('Execute the creation of the  GFF file of the genome FASTA.')
        self.pushbutton_execute.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_execute.clicked.connect(self.pushbutton_execute_clicked)

        # create and configure "pushbutton_close"
        pushbutton_close = QPushButton('Close')
        pushbutton_close.setToolTip('Cancel the creation of the  GFF file of the genome FASTA and close the window.')
        pushbutton_close.setCursor(QCursor(Qt.PointingHandCursor))
        pushbutton_close.clicked.connect(self.pushbutton_close_clicked)

        # create and configure "gridlayout_buttons"
        gridlayout_buttons = QGridLayout()
        gridlayout_buttons.setColumnStretch(0, 10)
        gridlayout_buttons.setColumnStretch(1, 1)
        gridlayout_buttons.setColumnStretch(2, 1)
        gridlayout_buttons.addWidget(self.pushbutton_execute, 0, 1, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(pushbutton_close, 0, 2, alignment=Qt.AlignCenter)

        # create and configure "groupbox_buttons"
        groupbox_buttons = QGroupBox()
        groupbox_buttons.setObjectName('groupbox_buttons')
        groupbox_buttons.setStyleSheet('QGroupBox#groupbox_buttons {border: 0px;}')
        groupbox_buttons.setLayout(gridlayout_buttons)

        # create and configure "gridlayout_central"
        gridlayout_central = QGridLayout()
        gridlayout_central.setRowStretch(0, 1)
        gridlayout_central.setRowStretch(1, 1)
        gridlayout_central.setRowStretch(2, 1)
        gridlayout_central.setRowStretch(3, 1)
        gridlayout_central.setRowStretch(4, 1)
        gridlayout_central.setColumnStretch(0, 0)
        gridlayout_central.setColumnStretch(1, 1)
        gridlayout_central.setColumnStretch(2, 0)
        gridlayout_central.addWidget(label_head, 0, 1)
        gridlayout_central.addWidget(QLabel(), 1, 1)
        gridlayout_central.addWidget(groupbox_data, 2, 1)
        gridlayout_central.addWidget(QLabel(), 3, 1)
        gridlayout_central.addWidget(groupbox_buttons, 4, 1)

        # create and configure "groupbox_central"
        groupbox_central = QGroupBox()
        groupbox_central.setLayout(gridlayout_central)

        # create and configure "vboxlayout"
        vboxlayout = QVBoxLayout(self)
        vboxlayout.addWidget(groupbox_central)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # set initial value in "lineedit_threads"
        self.lineedit_threads.setText('4')

        # populate data in "combobox_reference_species_name"
        self.combobox_reference_species_name_populate()

        # set initial value in "lineedit_target_genome_path"
        self.lineedit_target_genome_path.setText('')

    #---------------

    def check_inputs(self):
        '''
        Check the content of each input and do the actions linked to its value.
        '''

        # initialize the control variable
        OK = True

        # check "lineedit_threads" when the editing finished
        if not self.lineedit_threads_editing_finished():
            OK = False

        # check "lineedit_target_genome_path" when the editing finished
        if not self.lineedit_target_genome_path_editing_finished():
            OK = False

        # check all inputs are OK
        if OK:
            self.parent.statusBar().showMessage('')
        else:
            self.parent.statusBar().showMessage('WARNING: One or more input values are wrong or empty.')

        # enable "pushbutton_execute"
        if OK and self.lineedit_threads.text() != '' and self.lineedit_target_genome_path.text() != '':
            self.pushbutton_execute.setEnabled(True)
        else:
            self.pushbutton_execute.setEnabled(False)

        # return the control variable
        return OK

    #---------------

    def lineedit_threads_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_threads"
        '''

        # initialize the control variable
        OK = True

        # chek if "lineedit_threads" is empty
        if self.lineedit_threads.text() == '':
            OK = False
            self.lineedit_threads.setStyleSheet('background-color: white')

        # chek if "lineedit_threads" is an integer number between 1 and the CPUs number in the computer
        elif self.lineedit_threads.text() != '' and not genlib.check_int(self.lineedit_threads.text(), minimum=1, maximum=os.cpu_count()):
            OK = False
            self.lineedit_threads.setStyleSheet('background-color: red')
            text = f'The value of threads number has to be an integer number between 1 and {os.cpu_count()} (threads available in the computer).'
            QMessageBox.critical(self, self.title, text, buttons=QMessageBox.Ok)

        else:
            self.lineedit_threads.setStyleSheet('background-color: white')

        # return the control variable
        return OK

    #---------------

    def combobox_reference_species_name_populate(self):
        '''
        Populate data in "combobox_reference_species_name".
        '''

        # set the species names list
        species_names_list = [genlib.get_quercus_acutissima_name(), genlib.get_quercus_dentata_name(), genlib.get_quercus_glauca_name(), genlib.get_quercus_ilex_name(), genlib.get_quercus_lobata_name(), genlib.get_quercus_robur_name(), genlib.get_quercus_rubra_name(), genlib.get_quercus_suber_name(), genlib.get_quercus_variabilis_name()]

        # populate data in "combobox_reference_species_name"
        self.combobox_reference_species_name.addItems(sorted(species_names_list))

        # select the item corresponding to "Quercus lobata"
        index = self.combobox_reference_species_name.findText(genlib.get_quercus_lobata_name())
        if index != -1:
            self.combobox_reference_species_name.setCurrentIndex(index)

        # simultate "combobox_reference_species_name" index has changed
        self.combobox_reference_species_name_currentIndexChanged()

    #---------------

    def combobox_reference_species_name_currentIndexChanged(self):
        '''
        Process the event when an item of "combobox_reference_species_name" has been selected.
        '''

        # check the content of inputs
        self.check_inputs()

    #---------------

    def lineedit_target_genome_path_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_target_genome_path"
        '''

        # initialize the control variable
        OK = True

        # chek if "lineedit_target_genome_path" is empty
        if self.lineedit_target_genome_path.text() == '':
            self.lineedit_target_genome_path.setStyleSheet('background-color: white')
            OK = False

        # chek if "lineedit_target_genome_path" exists
        elif not os.path.isfile(self.lineedit_target_genome_path.text()):
            self.lineedit_target_genome_path.setStyleSheet('background-color: red')
            OK = False

        else:
            self.lineedit_target_genome_path.setStyleSheet('background-color: white')

        # return the control variable
        return OK

    #---------------

    def pushbutton_search_target_genome_path_clicked(self):
        '''
        Search and select the genome FASTA path.
        '''

        # search the genome FASTA patj
        (target_genome_path, _) = QFileDialog.getOpenFileName(self, f'{self.head} - Selection of genome FASTA file', os.path.expanduser('~'), "FASTA files (*.fasta *.FASTA *.fas *.FAS *.fa *.FA *.fsa *.FSA *.fna *.FNA *.faa *.FAA);;all (*.*)")

        # set "lineedit_target_genome_path" with the genome FASTA path selected
        if target_genome_path != '':
            self.lineedit_target_genome_path.setText(target_genome_path)

        # check the content of inputs
        self.check_inputs()

    #---------------

    def pushbutton_execute_clicked(self):
        '''
        Execute the process.
        '''

        # initialize the control variable
        OK = True

        # confirm the process is executed
        text = 'The process to create the genome GFF file is going to be run.\n\nAre you sure to continue?'
        botton = QMessageBox.question(self, self.title, text, buttons=QMessageBox.Yes|QMessageBox.No, defaultButton=QMessageBox.No)
        if botton == QMessageBox.No:
            OK = False

        # execute the process
        if OK:

            # get the threads number
            threads = self.lineedit_threads.text()

            # get the FASTA type
            fasta_type = self.combobox_reference_species_name.currentText()

            # get the genome FASTA path
            target_genome_path = self.lineedit_target_genome_path.text()
            if sys.platform.startswith('win32'):
                target_genome_path = genlib.windows_path_2_wsl_path(target_genome_path)

            # create and execute "DialogProcess"
            process = dialogs.DialogProcess(self, self.head, self.create_gff_file, threads, fasta_type, target_genome_path)
            process.exec()

        # close the windows
        # -- if OK:
        # --     self.pushbutton_close_clicked()

    #---------------

    def pushbutton_close_clicked(self):
        '''
        Close the window.
        '''

        self.parent.current_subwindow = None
        self.close()
        self.parent.set_background_image()

   #---------------

    def create_gff_file(self, process, threads, reference_species_name, target_genome_path):
        '''
        Run a process to create a genome GFF file.
        '''

        # initialize the control variable
        OK = True

        # warn that the log window does not have to be closed
        process.write('Do not close this window, please wait!\n')

        # warn that the requirements are being verified
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Checking process requirements ...\n')

        # check the quercusTOA config file
        if OK:
            if not os.path.isfile(genlib.get_app_config_file()):
                process.write(f'*** ERROR: The {genlib.get_app_short_name()} config file does not exist. Please, recreate it.\n')
                OK = False

        # warn that the requirements are OK
        if OK:
            process.write('Process requirements are OK.\n')

        # determine the temporal directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Determining the temporal directory ...\n')
            temp_dir = genlib.get_temp_dir()
            command = f'mkdir -p {temp_dir}'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write(f'The directory path is {temp_dir}.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # determine the run directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Determining the run directory ...\n')
            result_dir = self.app_config_dict['Environment parameters']['result_dir']
            current_run_dir = genlib.get_current_run_dir(result_dir, genlib.get_result_run_subdir(), genlib.get_process_create_gff_file_code())
            command = f'mkdir -p {current_run_dir}'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write(f'The directory path is {current_run_dir}.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # build the script
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            script_name = f'{genlib.get_process_create_gff_file_code()}-process.sh'
            process.write(f'Building the process script {script_name} ...\n')
            (OK, _) = self.build_create_gff_file_script(temp_dir, script_name, current_run_dir, threads, reference_species_name, target_genome_path)
            if OK:
                process.write('The file is built.\n')
            else:
                process.write('*** ERROR: The file could not be built.\n')

        # copy the script to the current run directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Copying the script {script_name} to the directory {current_run_dir} ...\n')
            command = f'cp {temp_dir}/{script_name} {current_run_dir}; [ $? -eq 0 ] &&  exit 0 || exit 1'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write('The file is copied.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')

        # set run permision to the script
        if OK and not sys.platform.startswith('win32'):
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Setting on the run permision of {script_name} ...\n')
            command = f'chmod u+x {current_run_dir}/{script_name}'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write('The run permision is set.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # build the starter script in the temporal directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            starter_name = f'{genlib.get_process_create_gff_file_code()}-process-starter.sh'
            process.write(f'Building the starter script {starter_name} ...\n')
            (OK, _) = genlib.build_starter(temp_dir, starter_name, script_name, current_run_dir)
            if OK:
                process.write('The file is built.\n')
            if not OK:
                process.write('***ERROR: The file could not be built.\n')

        # copy the starter script to the current run directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Copying the starter {starter_name} to the directory {current_run_dir} ...\n')
            command = f'cp {temp_dir}/{starter_name} {current_run_dir}; [ $? -eq 0 ] &&  exit 0 || exit 1'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write('The file is copied.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')

        # set run permision to the starter script
        if OK and not sys.platform.startswith('win32'):
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Setting on the run permision of {starter_name} ...\n')
            command = f'chmod u+x {current_run_dir}/{starter_name}'
            rc = genlib.run_command(command, process, is_script=False)
            if rc == 0:
                process.write('The run permision is set.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # submit the starter
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Submitting the starter {starter_name} ...\n')
            command = f'{current_run_dir}/{starter_name} &'
            rc = genlib.run_command(command, process, is_script=True)
            if rc == 0:
                process.write('The script is submitted.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # warn that the log window can be closed
        process.write(f'{genlib.get_separator()}\n')
        process.write('You can close this window now.\n')

        # return the control variable
        return OK

    #---------------

    def build_create_gff_file_script(self, directory, script_name, current_run_dir, threads, reference_species_name, target_genome_path):
        '''
        Build the script to create a GFF file.
        '''

        # initialize the control variable and error list
        OK = True
        error_list = []

        # get the Miniforge3 directory and its bin subdirectory
        miniforge3_dir = ''
        if sys.platform.startswith('win32'):
            miniforge3_dir = genlib.get_miniforge3_dir_in_wsl()
        elif sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
            miniforge3_dir = genlib.get_miniforge3_current_dir()
        miniforge3_bin_dir = f'{miniforge3_dir}/bin'

        # get items from dictionary of application configuration
        reference_genome_path = ''
        reference_gff_path = ''
        if reference_species_name == genlib.get_quercus_acutissima_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qacutissima_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qacutissima_gff_path']
        elif reference_species_name == genlib.get_quercus_dentata_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qdentata_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qdentata_gff_path']
        elif reference_species_name == genlib.get_quercus_glauca_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qglauca_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qglauca_gff_path']
        elif reference_species_name == genlib.get_quercus_ilex_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qilex_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qilex_gff_path']
        elif reference_species_name == genlib.get_quercus_lobata_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qlobata_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qlobata_gff_path']
        elif reference_species_name == genlib.get_quercus_robur_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qrobur_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qrobur_gff_path']
        elif reference_species_name == genlib.get_quercus_rubra_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qrubra_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qrubra_gff_path']
        elif reference_species_name == genlib.get_quercus_suber_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qsuber_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qsuber_gff_path']
        elif reference_species_name == genlib.get_quercus_variabilis_name():
            reference_genome_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qvariabilis_genome_path']
            reference_gff_path = self.app_config_dict[f'{genlib.get_app_short_name()} database']['qvariabilis_gff_path']

        # set the parameters file
        params_file = f'{current_run_dir}/{genlib.get_params_file_name()}'

        # set the temporal directories
        temp_dir = f'{current_run_dir}/temp'
        temp_liftoff_dir = f'{temp_dir}/temp-liftoff'

        #
        target_gff_file = f'{current_run_dir}/target.gff'
        unmapped_features_file = f'{temp_dir}/unmapped-features.txt'

        # set the script path
        script_path = f'{directory}/{script_name}'

        # write the script
        try:
            with open(script_path, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
                file_id.write( '#!/bin/bash\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write(f'export PATH={miniforge3_bin_dir}:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin\n')
                file_id.write( 'SEP="#########################################"\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write(f'STATUS_DIR={genlib.get_status_dir(current_run_dir)}\n')
                file_id.write(f'SCRIPT_STATUS_OK={genlib.get_status_ok(current_run_dir)}\n')
                file_id.write(f'SCRIPT_STATUS_WRONG={genlib.get_status_wrong(current_run_dir)}\n')
                file_id.write( 'mkdir -p $STATUS_DIR\n')
                file_id.write(f'mkdir -p {temp_dir}\n')
                file_id.write(f'mkdir -p {temp_liftoff_dir}\n')
                file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function init\n')
                file_id.write( '{\n')
                file_id.write( '    INIT_DATETIME=`date +%s`\n')
                file_id.write( '    FORMATTED_INIT_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME."\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write('function save_params\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Saving parameters ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/save-parms.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                file_id.write(f'        echo "[GFF creation parameters]" > {params_file}\n')
                file_id.write(f'        echo "target_genome_path = {target_genome_path}" >> {params_file}\n')
                file_id.write( '        RC=$?\n')
                file_id.write( '        if [ $RC -ne 0 ]; then manage_error echo $RC; fi\n')
                file_id.write( '        echo "Parameters are saved."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function create_gff_file\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Creating the GFF file ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/create-gff-file.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_liftoff_environment()}\n')
                file_id.write( '        /usr/bin/time \\\n')
                file_id.write( '            liftoff \\\n')
                file_id.write(f'                -p {threads} \\\n')
                file_id.write(f'                -g {reference_gff_path} \\\n')
                file_id.write(f'                -o {target_gff_file} \\\n')
                file_id.write(f'                -u {unmapped_features_file} \\\n')
                file_id.write(f'                -dir {temp_liftoff_dir} \\\n')
                file_id.write(f'                {target_genome_path} \\\n')
                file_id.write(f'                {reference_genome_path}\n')
                file_id.write( '        RC=$?\n')
                file_id.write( '        if [ $RC -ne 0 ]; then manage_error blastn $RC; fi\n')
                file_id.write( '        conda deactivate\n')
                file_id.write( '        echo "File is created."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function calculate_sequence_identity\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Calculating the sequence identity ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/calculate-sequence-identity.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_liftofftools_environment()}\n')
                file_id.write( '        /usr/bin/time \\\n')
                file_id.write( '            liftofftools \\\n')
                file_id.write( '                variants \\\n')
                file_id.write(f'                -r {reference_genome_path} \\\n')
                file_id.write(f'                -rg {reference_gff_path} \\\n')
                file_id.write(f'                -t {target_genome_path} \\\n')
                file_id.write(f'                -tg {target_gff_file}\n')
                file_id.write( '        RC=$?\n')
                file_id.write( '        if [ $RC -ne 0 ]; then manage_error blastn $RC; fi\n')
                file_id.write( '        conda deactivate\n')
                file_id.write( '        echo "Calculation is done."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function compare_gene_order\n')
                file_id.write( '{\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Comparing the gene order ..."\n')
                file_id.write(f'    cd {current_run_dir}\n')
                file_id.write( '    STEP_STATUS=$STATUS_DIR/compare-gene-order.ok\n')
                file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                file_id.write( '        echo "This step was previously run."\n')
                file_id.write( '    else\n')
                file_id.write(f'        source {miniforge3_bin_dir}/activate {genlib.get_liftofftools_environment()}\n')
                file_id.write( '        /usr/bin/time \\\n')
                file_id.write( '            liftofftools \\\n')
                file_id.write( '                synteny \\\n')
                file_id.write(f'                -r {reference_genome_path} \\\n')
                file_id.write(f'                -rg {reference_gff_path} \\\n')
                file_id.write(f'                -t {target_genome_path} \\\n')
                file_id.write(f'                -tg {target_gff_file}\n')
                file_id.write( '        RC=$?\n')
                file_id.write( '        if [ $RC -ne 0 ]; then manage_error blastn $RC; fi\n')
                file_id.write( '        conda deactivate\n')
                file_id.write( '        echo "Calculation is done."\n')
                file_id.write( '        touch $STEP_STATUS\n')
                file_id.write( '    fi\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function end\n')
                file_id.write( '{\n')
                file_id.write( '    END_DATETIME=`date +%s`\n')
                file_id.write( '    FORMATTED_END_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`\n')
                file_id.write( '    calculate_duration\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "Script ended OK at $FORMATTED_END_DATETIME with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    touch $SCRIPT_STATUS_OK\n')
                file_id.write( '    exit 0\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function manage_error\n')
                file_id.write( '{\n')
                file_id.write( '    END_DATETIME=`date +%s`\n')
                file_id.write( '    FORMATTED_END_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`\n')
                file_id.write( '    calculate_duration\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    echo "ERROR: $1 returned error $2"\n')
                file_id.write( '    echo "Script ended WRONG at $FORMATTED_END_DATETIME with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
                file_id.write( '    echo "$SEP"\n')
                file_id.write( '    touch $SCRIPT_STATUS_WRONG\n')
                file_id.write( '    exit 3\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'function calculate_duration\n')
                file_id.write( '{\n')
                file_id.write( '    DURATION=`expr $END_DATETIME - $INIT_DATETIME`\n')
                file_id.write( '    HH=`expr $DURATION / 3600`\n')
                file_id.write( '    MM=`expr $DURATION % 3600 / 60`\n')
                file_id.write( '    SS=`expr $DURATION % 60`\n')
                file_id.write( '    FORMATTED_DURATION=`printf "%03d:%02d:%02d\\n" $HH $MM $SS`\n')
                file_id.write( '}\n')
                file_id.write( '#-------------------------------------------------------------------------------\n')
                file_id.write( 'init\n')
                file_id.write( 'save_params\n')
                file_id.write( 'create_gff_file\n')
                file_id.write( 'calculate_sequence_identity\n')
                file_id.write( 'compare_gene_order\n')
                file_id.write( 'end\n')
        except Exception as e:
            error_list.append(f'*** EXCEPTION: "{e}".')
            error_list.append(f'*** ERROR: The file {script_path} is not created.')
            OK = False

        # return the control variable and error list
        return (OK, error_list)

    #---------------

#-------------------------------------------------------------------------------

class FormRestartGffCreation(QWidget):
    '''
    Class used to restart a GFF file creation.
    '''

    #---------------

    def __init__(self, parent):
        '''
        Create a class instance.
        '''

        # save parameters in instance variables
        self.parent = parent

        # call the init method of the parent class
        super().__init__()

        # set the dimensions window
        self.window_height = self.parent.WINDOW_HEIGHT - 100
        self.window_width = self.parent.WINDOW_WIDTH - 50

        # set the head and title
        self.head = 'Restart a GFF file creation'
        self.title = f'{genlib.get_app_short_name()} - {self.head}'

        # get the dictionary of application configuration
        self.app_config_dict = genlib.get_config_dict(genlib.get_app_config_file())

        # build the graphic user interface of the window
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

        # show the window
        self.show()

    #---------------

    def build_gui(self):
        '''
        Build the graphic user interface of the window.
        '''

        # set the width and height of the window
        self.setFixedSize(self.window_width, self.window_height)

        # move the window at center
        rectangle = self.frameGeometry()
        central_point = QGuiApplication.primaryScreen().availableGeometry().center()
        rectangle.moveCenter(central_point)
        self.move(rectangle.topLeft())

        # get font metrics information
        fontmetrics = QFontMetrics(QApplication.font())

        # create and configure "label_head"
        label_head = QLabel(self.head, alignment=Qt.AlignCenter)
        label_head.setStyleSheet('font: bold 14px; color: black; background-color: lightGray; max-height: 30px')

        # create and configure "tablewidget"
        self.tablewidget = QTableWidget()
        self.tablewidget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tablewidget.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.column_name_list = ['Process', 'Result dataset', 'Date', 'Time', 'Status']
        self.tablewidget.setColumnCount(len(self.column_name_list))
        self.tablewidget.setHorizontalHeaderLabels(self.column_name_list)
        self.tablewidget.setColumnWidth(0, 230)
        self.tablewidget.setColumnWidth(1, 280)
        self.tablewidget.setColumnWidth(2, 85)
        self.tablewidget.setColumnWidth(3, 70)
        self.tablewidget.setColumnWidth(4, 90)
        self.tablewidget.verticalHeader().setVisible(True)
        self.tablewidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.tablewidget.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tablewidget.currentCellChanged.connect(self.tablewidget_currentCellChanged)
        self.tablewidget.cellClicked.connect(self.tablewidget_cellClicked)
        self.tablewidget.cellDoubleClicked.connect(self.tablewidget_cellDoubleClicked)

        # create and configure "label_target_genome_path"
        label_target_genome_path = QLabel()
        label_target_genome_path.setText('Genome FASTA file')
        label_target_genome_path.setFixedWidth(fontmetrics.width('9'*20))

        # create and configure "lineedit_target_genome_path"
        self.lineedit_target_genome_path = QLineEdit()
        self.lineedit_target_genome_path.editingFinished.connect(self.check_inputs)
        self.lineedit_target_genome_path.setDisabled(True)

        # create and configure "gridlayout_data"
        gridlayout_data = QGridLayout()
        gridlayout_data.setColumnStretch(0, 1)
        gridlayout_data.setColumnStretch(1, 1)
        gridlayout_data.addWidget(self.tablewidget, 0, 0, 1, 11)
        gridlayout_data.addWidget(label_target_genome_path, 1, 0, 1, 1)
        gridlayout_data.addWidget(self.lineedit_target_genome_path, 1, 1, 2, 7)

        # create and configure "groupbox_data"
        groupbox_data = QGroupBox()
        groupbox_data.setObjectName('groupbox_data')
        groupbox_data.setStyleSheet('QGroupBox#groupbox_data {border: 0px;}')
        groupbox_data.setLayout(gridlayout_data)

        # create and configure "pushbutton_refresh"
        self.pushbutton_refresh = QPushButton('Refresh')
        self.pushbutton_refresh.setToolTip('Update the process list.')
        self.pushbutton_refresh.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_refresh.clicked.connect(self.pushbutton_refresh_clicked)

        # create and configure "pushbutton_execute"
        self.pushbutton_execute = QPushButton('Execute')
        self.pushbutton_execute.setToolTip('Restart the process selected.')
        self.pushbutton_execute.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushbutton_execute.clicked.connect(self.pushbutton_execute_clicked)

        # create and configure "pushbutton_close"
        pushbutton_close = QPushButton('Close')
        pushbutton_close.setToolTip('Close the window.')
        pushbutton_close.setCursor(QCursor(Qt.PointingHandCursor))
        pushbutton_close.clicked.connect(self.pushbutton_close_clicked)

        # create and configure "gridlayout_buttons"
        gridlayout_buttons = QGridLayout()
        gridlayout_buttons.setColumnStretch(0, 15)
        gridlayout_buttons.setColumnStretch(1, 1)
        gridlayout_buttons.setColumnStretch(2, 1)
        gridlayout_buttons.setColumnStretch(3, 1)
        gridlayout_buttons.addWidget(self.pushbutton_refresh, 0, 1, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(self.pushbutton_execute, 0, 2, alignment=Qt.AlignCenter)
        gridlayout_buttons.addWidget(pushbutton_close, 0, 3, alignment=Qt.AlignCenter)

        # create and configure "groupbox_buttons"
        groupbox_buttons = QGroupBox()
        groupbox_buttons.setObjectName('groupbox_buttons')
        groupbox_buttons.setStyleSheet('QGroupBox#groupbox_buttons {border: 0px;}')
        groupbox_buttons.setLayout(gridlayout_buttons)

        # create and configure "gridlayout_central"
        gridlayout_central = QGridLayout()
        gridlayout_central.setRowStretch(0, 1)
        gridlayout_central.setRowStretch(1, 1)
        gridlayout_central.setRowStretch(2, 10)
        gridlayout_central.setRowStretch(3, 1)
        gridlayout_central.setColumnStretch(0, 1)
        gridlayout_central.addWidget(label_head, 0, 0)
        gridlayout_central.addWidget(QLabel(), 1, 0)
        gridlayout_central.addWidget(groupbox_data, 2, 0)
        gridlayout_central.addWidget(groupbox_buttons, 3, 0)

        # create and configure "groupbox_central"
        groupbox_central = QGroupBox()
        groupbox_central.setLayout(gridlayout_central)

        # create and configure "vboxlayout"
        vboxlayout = QVBoxLayout(self)
        vboxlayout.addWidget(groupbox_central)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load data in "tablewidget"
        self.load_tablewidget()

        # initialize "lineedit_ftarget_genome_path"
        self.lineedit_target_genome_path.setText('')

    #---------------

    def check_inputs(self):
        '''
        Check the content of each input and do the actions linked to its value.
        '''

        # initialize the control variable
        OK = True

        # get the list of rows selected
        row_list = []
        for idx in self.tablewidget.selectionModel().selectedIndexes():
            row_list.append(idx.row())
        row_list = list(set(row_list))

        # check "lineedit_target_genome_path" when the editing finished
        if not self.lineedit_target_genome_path_editing_finished():
            OK = False

        # check all inputs are OK
        if OK:
            self.parent.statusBar().showMessage('')
        else:
            self.parent.statusBar().showMessage('There are one or more inputs without data or with wrong value.')

        # enable "pushbutton_execute"
        if OK and len(row_list) == 1 and self.lineedit_target_genome_path.text() != '':
            self.pushbutton_execute.setEnabled(True)
        else:
            self.pushbutton_execute.setEnabled(False)
            OK = False

        # return the control variable
        return OK

    #---------------

    def tablewidget_currentCellChanged(self, row, _):
        '''
        Perform necessary actions after changing the current "tablewidget" cell.
        '''

        # check if there is a row selected
        if row >= 0:

            # get the dictionary of application configuration
            result_dir = self.app_config_dict['Environment parameters']['result_dir']

            # get the result dataset corresponding to the row clicked
            result_dataset = self.tablewidget.item(row, 1).text()

            # get the parameters file corresponding to the result dataset
            params_file = f'{result_dir}{os.sep}{genlib.get_result_run_subdir()}{os.sep}{result_dataset}{os.sep}{genlib.get_params_file_name()}'
            if sys.platform.startswith('win32'):
                params_file = genlib.wsl_path_2_windows_path(params_file)

            # get the dictionary of parameters
            params_dict = genlib.get_config_dict(params_file)

            # set the genome FASTA path in "lineedit_target_genome_path"
            target_genome_path = params_dict['GFF creation parameters']['target_genome_path']
            if sys.platform.startswith('win32'):
                target_genome_path = genlib.wsl_path_2_windows_path(target_genome_path)
            self.lineedit_target_genome_path.setText(target_genome_path)

        # check the content of inputs
        self.check_inputs()

    #---------------

    def tablewidget_cellClicked(self, _, __):
        '''
        Perform necessary actions after clicking on a "tablewidget" cell.
        '''

        # check the content of inputs
        self.check_inputs()

    #---------------

    def tablewidget_cellDoubleClicked(self, _, __):
        '''
        Perform necessary actions after double clicking on "tablewidget" cell.
        '''

        # check the content of inputs
        OK = self.check_inputs()

        # if inputs are OK, simulate a click on "pushbutton_execute"
        if OK:
            self.pushbutton_execute_clicked()

    #---------------

    def lineedit_target_genome_path_editing_finished(self):
        '''
        Perform necessary actions after finishing editing "lineedit_target_genome_path"
        '''

        # initialize the control variable
        OK = True

        # return the control variable
        return OK

    #---------------

    def pushbutton_refresh_clicked(self):
        '''
        Refresh "tablewidget" and initialize other inputs.
        '''

        # check the content of inputs
        self.initialize_inputs()

        # check the content of inputs
        self.check_inputs()

    #---------------

    def pushbutton_execute_clicked(self):
        '''
        Execute the process.
        '''

        # initialize the control variable
        OK = True

        # get the list of rows selected
        row_list = []
        for idx in self.tablewidget.selectionModel().selectedIndexes():
            row_list.append(idx.row())
        row_list = list(set(row_list))

        # check if there is only a row selected
        if len(row_list) != 1:
            title = f'{genlib.get_app_short_name()} - {self.head}'
            text = 'One row has to be selected.'
            QMessageBox.critical(self, title, text, buttons=QMessageBox.Ok)
            OK = False

        # confirm the process is executed
        if OK:
            text = 'The GFF file creation is going to be restarted.\n\nAre you sure to continue?'
            botton = QMessageBox.question(self, self.title, text, buttons=QMessageBox.Yes|QMessageBox.No, defaultButton=QMessageBox.No)
            if botton == QMessageBox.No:
                OK = False

        # execute the process
        if OK:

            # get the result dataset
            result_dataset = self.tablewidget.item(row_list[0], 1).text()

            # create and execute "DialogProcess"
            process = dialogs.DialogProcess(self, self.head, self.restart_gff_file_creation, result_dataset)
            process.exec()

        # close the windows
        # -- if OK:
        # --     self.pushbutton_close_clicked()

    #---------------

    def pushbutton_close_clicked(self):
        '''
        Close the window.
        '''

        self.parent.current_subwindow = None
        self.close()
        self.parent.set_background_image()

    #---------------

    def load_tablewidget(self):
        '''
        Load data in "tablewidget".
        '''

        # get the result directory
        result_dir = self.app_config_dict['Environment parameters']['result_dir']

        # set the type, name and code of the annotation pipeline datasets
        process_type = genlib.get_result_run_subdir()
        process_name = genlib.get_process_create_gff_file_name()
        process_code = genlib.get_process_id(process_name)

        # get the process dictionary
        process_dict = genlib.get_process_dict()

        # get the log directory
        log_dir = f'{result_dir}/{process_type}'
        if sys.platform.startswith('win32'):
            log_dir = genlib.wsl_path_2_windows_path(log_dir)

        # set the command to get the result datasets of annotation pipeline in the log directory
        command = ''
        if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
            if process_name == 'all':
                command = f'ls -d {log_dir}/*  | xargs -n 1 basename'
            else:
                command = f'ls -d {log_dir}/{process_code}-*  | xargs -n 1 basename'
        elif sys.platform.startswith('win32'):
            log_dir = log_dir.replace('/', '\\')
            if process_name == 'all':
                command = f'dir /a:d /b {log_dir}\\*'
            else:
                command = f'dir /a:d /b {log_dir}\\{process_code}-*'

        # run the command to get the result datasets of enrichment analysis in the log directory
        output = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=False)

        # initialize the result dataset dictionary
        result_dataset_dict = {}

        # build the result dataset dictionary
        for line in output.stdout.split('\n'):
            if line != '':


                # get data
                result_dataset_id = line.strip()
                try:
                    pattern = r'^(.+)\-(.+)\-(.+)$'
                    mo = re.search(pattern, result_dataset_id)
                    process_code = mo.group(1).strip()
                    process_name = process_dict[process_code]['name']
                    yymmdd = mo.group(2)
                    hhmmss = mo.group(3)
                    date = f'20{yymmdd[:2]}-{yymmdd[2:4]}-{yymmdd[4:]}'
                    time = f'{hhmmss[:2]}:{hhmmss[2:4]}:{hhmmss[4:]}'
                except:    # pylint: disable=bare-except
                    process_name = 'unknown process'
                    date = '0000-00-00'
                    time = '00:00:00'

                # determine the status
                status_ok = os.path.isfile(genlib.get_status_ok(os.path.join(log_dir, result_dataset_id)))
                status_wrong = os.path.isfile(genlib.get_status_wrong(os.path.join(log_dir, result_dataset_id)))
                status = ''
                if status_ok and not status_wrong:
                    status = 'OK'
                elif not status_ok and status_wrong:
                    status = 'wrong'
                elif not status_ok and not status_wrong:
                    status = 'not finished'
                elif status_ok and status_wrong:
                    status = 'undetermined'

                # insert data in the dictionary when the the status is wrong
                if status == 'wrong':
                    key = f'{process_name}-{result_dataset_id}'
                    result_dataset_dict[key] = {'process': process_name, 'result_dataset_id': result_dataset_id, 'date': date, 'time': time, 'status': status}

        # initialize "tablewidget"
        self.tablewidget.clearContents()

        # set the rows number of "tablewidget"
        self.tablewidget.setRowCount(len(result_dataset_dict))

        # load data in "tablewidget" for the OK result datasets of annotation pipeline
        if not result_dataset_dict:
            text = 'There is no run ended wrong.'
            QMessageBox.warning(self, self.title, text, buttons=QMessageBox.Ok)
        else:
            row = 0
            for key in sorted(result_dataset_dict.keys()):
                self.tablewidget.setItem(row, 0, QTableWidgetItem(result_dataset_dict[key]['process']))
                self.tablewidget.setItem(row, 1, QTableWidgetItem(result_dataset_dict[key]['result_dataset_id']))
                self.tablewidget.setItem(row, 2, QTableWidgetItem(result_dataset_dict[key]['date']))
                self.tablewidget.setItem(row, 3, QTableWidgetItem(result_dataset_dict[key]['time']))
                self.tablewidget.setItem(row, 4, QTableWidgetItem(result_dataset_dict[key]['status']))
                row += 1

   #---------------

    def restart_gff_file_creation(self, process, result_dataset):
        '''
        Restart a GFF file creation.
        '''

        # initialize the control variable
        OK = True

        # warn that the log window does not have to be closed
        process.write('Do not close this window, please wait!\n')

        # warn that the requirements are being verified
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Checking process requirements ...\n')

        # check the quercusTOA config file
        if OK:
            if not os.path.isfile(genlib.get_app_config_file()):
                process.write(f'*** ERROR: The {genlib.get_app_short_name()} config file does not exist. Please, recreate it.\n')
                OK = False

        # warn that the requirements are OK
        if OK:
            process.write('Process requirements are OK.\n')

        # determine the run directory
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write('Determining the run directory ...\n')
            result_dir = self.app_config_dict['Environment parameters']['result_dir']
            current_run_dir = f'{result_dir}/{genlib.get_result_run_subdir()}/{result_dataset}'
            process.write(f'The directory path is {current_run_dir}.\n')

        # set the starter script name
        if OK:
            starter_name = f'{genlib.get_process_create_gff_file_code()}-process-starter.sh'

        # submit the starter
        if OK:
            process.write(f'{genlib.get_separator()}\n')
            process.write(f'Submitting the starter {starter_name} ...\n')
            command = f'{current_run_dir}/{starter_name} &'
            rc = genlib.run_command(command, process, is_script=True)
            if rc == 0:
                process.write('The script is submitted.\n')
            else:
                process.write(f'*** ERROR: RC {rc} in command -> {command}\n')
                OK = False

        # warn that the log window can be closed
        process.write(f'{genlib.get_separator()}\n')
        process.write('You can close this window now.\n')

        # return the control variable
        return OK

    #---------------

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print(f'This file contains the classes related to the comparative genomics used in {genlib.get_app_long_name()}')
    sys.exit(0)

#-------------------------------------------------------------------------------

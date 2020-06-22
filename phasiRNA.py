from PyQt5 import QtWidgets, QtCore, QtGui
from fuzzysearch import find_near_matches
import math
import collections
import string
import re
import sys


class Ui_phasiRNA(object):
    def setupUi(self, phasiRNA):
        phasiRNA.setObjectName("phasiRNA")
        phasiRNA.resize(540, 386)
        self.centralwidget = QtWidgets.QWidget(phasiRNA)
        self.centralwidget.setObjectName("centralwidget")
        self.file_button = QtWidgets.QPushButton(self.centralwidget)
        self.file_button.setGeometry(QtCore.QRect(50, 40, 151, 41))
        font = QtGui.QFont()
        font.setFamily("Agency FB")
        font.setPointSize(9)
        self.file_button.setFont(font)
        self.file_button.setObjectName("file_button")
        self.goalfile_button = QtWidgets.QPushButton(self.centralwidget)
        self.goalfile_button.setGeometry(QtCore.QRect(50, 110, 151, 41))
        font = QtGui.QFont()
        font.setFamily("Agency FB")
        font.setPointSize(9)
        self.goalfile_button.setFont(font)
        self.goalfile_button.setObjectName("goalfile_button")
        self.run_button = QtWidgets.QPushButton(self.centralwidget)
        self.run_button.setGeometry(QtCore.QRect(50, 270, 151, 41))
        font = QtGui.QFont()
        font.setFamily("Agency FB")
        font.setPointSize(9)
        self.run_button.setFont(font)
        self.run_button.setObjectName("run_button")
        self.file_input = QtWidgets.QTextBrowser(self.centralwidget)
        self.file_input.setGeometry(QtCore.QRect(230, 40, 281, 41))
        self.file_input.setObjectName("file_input")
        self.goalfile_input = QtWidgets.QTextBrowser(self.centralwidget)
        self.goalfile_input.setGeometry(QtCore.QRect(230, 110, 281, 41))
        self.goalfile_input.setObjectName("goalfile_input")
        self.run_output = QtWidgets.QTextBrowser(self.centralwidget)
        self.run_output.setGeometry(QtCore.QRect(230, 250, 281, 81))
        self.run_output.setObjectName("run_output")
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(50, 180, 151, 41))
        font = QtGui.QFont()
        font.setFamily("Agency FB")
        font.setPointSize(9)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setGeometry(QtCore.QRect(230, 180, 281, 41))
        self.textBrowser.setObjectName("textBrowser")
        phasiRNA.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(phasiRNA)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 540, 18))
        self.menubar.setObjectName("menubar")
        phasiRNA.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(phasiRNA)
        self.statusbar.setObjectName("statusbar")
        phasiRNA.setStatusBar(self.statusbar)

        self.retranslateUi(phasiRNA)
        QtCore.QMetaObject.connectSlotsByName(phasiRNA)

    def retranslateUi(self, phasiRNA):
        _translate = QtCore.QCoreApplication.translate
        phasiRNA.setWindowTitle(_translate("phasiRNA", "phasiRNA"))
        self.file_button.setText(_translate("phasiRNA", "选择文件"))
        self.goalfile_button.setText(_translate("phasiRNA", "选择参考序列文件"))
        self.run_button.setText(_translate("phasiRNA", "运行"))
        self.pushButton.setText(_translate("phasiRNA", "选择输出文件"))


class MainUi(QtWidgets.QMainWindow, Ui_phasiRNA):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_phasiRNA.__init__(self)
        self.setupUi(self)
        self.file_button.clicked.connect(self.selfile)
        self.goalfile_button.clicked.connect(self.selgoalfile)
        self.pushButton.clicked.connect(self.selfileout)
        self.run_button.clicked.connect(self.run)
        self.fileDialog = QtWidgets.QFileDialog(self)

    def selfile(self):
        self.selfile_input, _filter = self.fileDialog.getOpenFileName()
        self.file_input.setText("".join(self.selfile_input) + "\n")
        print("input" + "".join(self.selfile_input))

    def selgoalfile(self):
        self.selgoalfile_input, _filter = self.fileDialog.getOpenFileName()
        self.goalfile_input.setText("".join(self.selgoalfile_input))
        print("input" + "".join(self.selgoalfile_input))

    def selfileout(self):
        self.file_out, _filter = self.fileDialog.getSaveFileName()
        self.textBrowser.setText("".join(self.file_out))
        print("output" + "".join(self.file_out))

    def run(self):
        out = "Job Starts" + "\n"
        self.run_output.setText(out)
        PR = PhasiRna
        f2f = PR.fastq2fasta
        ref = PR.refseq
        rev = PR.reverseComplement
        any = PR.anyls
        result = PR.main(f2f, self.selfile_input, self.selgoalfile_input, ref, any, rev)
        output = PR.writeoutput(self.selgoalfile_input, result, self.file_out)
        out += "All Job Done" + "\n"
        self.run_output.setText(out)
        print("All Job Done")


class PhasiRna:

    def fastq2fasta(input_file):
        reads = []
        F = open(input_file, 'r')
        line_count = 0
        for line in F:
            line_count += 1
            if line_count % 4 == 2:
                reads.append(line.rstrip())
        return reads

    def refseq(input_file):
        seq = []
        sequences = open(input_file, 'r')
        for line in sequences:
            if line.startswith(">"):
                line1 = sequences.readline().strip()
                seq.append(line1)
        return seq

    def reverseComplement(sequence):
        complement = str.maketrans('ATCGN', 'TAGCN')
        return sequence.upper().translate(complement)[::-1]

    def anyls(a):
        dic = collections.Counter(a)
        return dic.most_common()

    def main(fastq2fasta, seqfile, reffile, refseq, anyls, reverseComplement):
        reads = fastq2fasta(seqfile)
        seq = refseq(reffile)
        reads_abundance = anyls(reads)
        score_list = []

        reads_dict = {}
        for i in reads_abundance:
            key = [i[0]]
            value = [i[1]]
            d = dict(zip(key, value))
            reads_dict.update(d)
        keys = list(reads_dict.keys())

        reversereads = []
        for i in keys:
            reverseq = reverseComplement(i)
            reversereads.append(reverseq)

        for b in seq:
            reads_align = {}
            for m in keys:
                c = find_near_matches(m, b, max_substitutions=2, max_deletions=0, max_insertions=0)
                sites = []
                for i in c:
                    sites.append(i.start + 1)
                reads_align[m] = sites

            reversereads_align = {}
            for n in reversereads:
                d = find_near_matches(n, b, max_substitutions=2, max_deletions=0, max_insertions=0)
                reversesites = []
                for i in d:
                    reversesites.append(i.start + 3)
                reversereads_align[n] = reversesites

            values1 = list(reads_align.values())
            values2 = list(reversereads_align.values())

            values = []
            if len(reads[0]) == 21:
                for (i, j) in zip(values1, values2):
                    z = list(set(i + j))
                    bin_sites = []
                    for k in z:
                        if k % 21 == 0:
                            bin_site = 21
                        else:
                            bin_site = k % 21
                        bin_sites.append(bin_site)
                    values.append(bin_sites)

                alignment = dict(zip(keys, values))
                alignment_list = list(alignment.items())
                align_length = len(alignment_list)

                new_list = []
                for i in range(align_length):
                    if not alignment_list[i][-1]:
                        pass
                    else:
                        new_list.append(alignment_list[i])

                number_list = []
                for i in range(1, 22):
                    number_list.append(i)

            elif len(reads[0]) == 24:
                for (i, j) in zip(values1, values2):
                    z = list(set(i + j))
                    bin_sites = []
                    for k in z:
                        if k % 24 == 0:
                            bin_site = 24
                        else:
                            bin_site = k % 24
                        bin_sites.append(bin_site)
                    values.append(bin_sites)

                alignment = dict(zip(keys, values))
                alignment_list = list(alignment.items())
                align_length = len(alignment_list)

                new_list = []
                for i in range(align_length):
                    if not alignment_list[i][-1]:
                        pass
                    else:
                        new_list.append(alignment_list[i])

                number_list = []
                for i in range(1, 25):
                    number_list.append(i)

            all_dict = {}
            for c in number_list:
                align_list = []
                for i in new_list:
                    try:
                        if c in i[1]:
                            align_list.append(i)
                    except:
                        if [c] == i[1]:
                            align_list.append(i)
                all_dict[c] = align_list

            phased_dict = {}
            for i in number_list:
                bin_abundance = 0
                for j in all_dict[i]:
                    bin_abundance += reads_dict[j[0]]
                phased_dict[i] = bin_abundance

            cluster_abundance = 0
            for i in phased_dict.values():
                cluster_abundance += i
            for key, value in phased_dict.items():
                if value == max(phased_dict.values()):
                    max_key = key
            phased_bin = list(phased_dict.keys())[list(phased_dict.values()).index(max(phased_dict.values()))]
            phased_ratio = max(phased_dict.values()) / cluster_abundance
            phased_number = 0
            for i in all_dict[phased_bin]:
                phased_number += 1
            phased_abundance = max(phased_dict.values())
            phased_score = phased_ratio * phased_number * math.log(phased_abundance)
            phased_score = round(phased_score, 3)
            score_list.append(phased_score)
        return score_list

    def writeoutput(file, score_li, outfile):
        fi = open(file, "r").readlines()
        name_list = []
        for line in fi:
            if line.startswith(">"):
                line = line.lstrip(">").rstrip()
                name_list.append(line)
        for i in range(len(score_li)):
            name_list[i] = name_list[i] + "\t" + str(score_li[i]) + "\n"
        file_out = open(outfile, "w")
        for c in name_list:
            file_out.write(c)
        file_out.close()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MainUi()
    window.show()
    sys.exit(app.exec_())

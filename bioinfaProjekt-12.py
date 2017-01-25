# -*- coding: utf-8 -*-

from cStringIO import StringIO
from Bio import Phylo
from Bio.Nexus import Trees
from ete2 import Tree
import os
import copy


class TreeData:
    def __init__(self, filename, data):
        self.filename = filename
        self.data = data
        self.phyloTree = Phylo.read(StringIO(data), "newick")
        self.eteTree = Tree(data + ';', format=1)


treeFiles = []


def menu():
    print("##############################################\n")
    print("Zadanie 12 - Kalkulator drzew filogenetycznych\n")
    print("Grupa:")
    print("Kamil Dąbrówka")
    print("Michał Kokocha")
    print("Mateusz Kowalkowski\n\n")
    print("--------MENU--------")
    print("1 - Wypisz drzewa")
    print("2 - Drzewo konsensusu")
    print("3 - Dystans RF")
    print("4 - Odcięcie")
    print("5 - Konwersja drzewo => rodzina zgodnych rozbić")
    print("\n0: Wyjście\n")
    print("##############################################\n\n")
    return


def readTree():
    print("WCZYTYWANIE DRZEW...")
    for file in os.listdir("trees"):
        if file.endswith(".newick"):
            f = open("trees/" + file, "r")
            treeFile = f.read()
            print(treeFile)
            handle = StringIO(treeFile)
            filename = os.path.splitext(file)[0]
            treeFiles.append(TreeData(filename, treeFile))
            f.close()
    print("DRZEWA WCZYTANE")
    return


def printTrees():
    print("Wypisywanie drzew\n")
    i = 0
    for tree in treeFiles:
        print("DRZEWO: " + str(i + 1))
        i = i + 1
        print tree.eteTree
    return;


def cutLeaf():
    print("Obcięcie podanego drzewa do drzewa filogenetycznego do zadanego podzbioru liści\n")

    tree_num = '-1'

    while (True):
        tree_num = raw_input("podaj drzewo do obciecia: ")
        try:
            tree_num = int(tree_num)
        except ValueError:
            print("Błąd: Niepoprawny numer drzewa lub znak. Spróbuj ponownie...")
            continue

        if (tree_num - 1 >= len(treeFiles) or tree_num < 1):
            print("Błąd: Nie ma takiego drzewa.")
        else:
            break

    print("Podane drzewo: ")

    print treeFiles[tree_num - 1].eteTree

    leafs = []
    for clade in treeFiles[tree_num - 1].phyloTree.find_clades():
        if clade.name is not None:
            leafs.append(str(clade.name))

    leafs = ''.join(leafs)
    leafs_given_not_valid = True

    while leafs_given_not_valid:
        leafs_not_to_cut = raw_input("Podaj zadany podzbiór lisci do ktorego obcięte zostanie drzewo: ")
        whole_string_flag = True
        for c in leafs_not_to_cut:
            if c not in leafs:
                whole_string_flag = False
                print("Błąd: Podany zbiór liści zawiera liście nieistniejące w wybranym drzewie")
                break

        if whole_string_flag:
            if len(leafs_not_to_cut) >= 2:
                leafs_given_not_valid = False
            else:
                print("Błąd: trzeba zostawić w drzewie co najmniej 2 liście")

    leafs_to_cut = []
    for leaf in leafs:
        if leaf not in leafs_not_to_cut:
            leafs_to_cut.append(leaf)

    leafs_to_cut = ''.join(leafs_to_cut)
    cutted_tree = treeFiles[tree_num - 1].phyloTree

    for i in leafs_to_cut:
        for clade in cutted_tree.find_clades():
            if str(clade.name) in leafs_to_cut and clade.is_terminal():
                cutted_tree.prune(clade)
    print("Wynikowe drzewo po obcięciu: ")
    Phylo.draw_ascii(cutted_tree)

    return


def consensus():
    print("Znajdywanie drzewa konsensusu dla wskazanych drzew\n")
    tree_num = '-1'
    con_trees = []
    while (True):
        tree_num = raw_input("Podaj drzewo lub X, jeśli skończyłeś wczytywanie: ")
        try:
            tree_num = int(tree_num)
        except ValueError:
            if tree_num == 'X':
                break
            else:
                print("Niepoprawny numer drzewa lub znak. Spróbuj ponownie...")
                continue

        if (tree_num - 1 >= len(treeFiles) or tree_num < 1):
            print("Błąd: Nie ma takiego drzewa.")
        else:
            if con_trees:  # jeśli już jakieś drzewo zostało wczytane
                if treeFiles[tree_num - 1].phyloTree.count_terminals() != con_trees[
                    0].phyloTree.count_terminals():  # sprawdź czy ma taką samą liczbę liści
                    print("Podane drzewo ma inną liczbę liści niż poprzednie")
                else:  # jeśli tak to dodaj
                    con_trees.append(treeFiles[tree_num - 1])
                    print("Dodano drzewo " + str(tree_num))

                    print treeFiles[tree_num - 1].eteTree

            else:  # jeśli to pierwsze wczytywane drzewo to dodaj
                con_trees.append(treeFiles[tree_num - 1])
                print("Dodane drzewo: " + str(tree_num))

                print treeFiles[tree_num - 1].eteTree

    threshold = -1.0
    while (threshold <= 0.0 or threshold > 1.0):
        threshold = float(raw_input("Podaj poziom toleracji dla drzewa konsensusu z zakresu (0, 1>: "))

    trees_to_con = []
    for t in treeFiles:
        tr = Trees.Tree(t.data)
        trees_to_con.append(tr)

    consensus_tree = Trees.consensus(trees_to_con, threshold)
    handle = StringIO(consensus_tree.to_string(plain_newick=True))
    tree = Phylo.read(handle, 'newick')
    Phylo.draw_ascii(tree)


def distance():
    print("Liczenie dystansu topologicznego RF pomiędzy parą drzew\n")
    tree_num = '-1'
    t1, t2 = None, None
    while (t1 == None or t2 == None):
        if not t1:
            tree_num = raw_input("Podaj pierwsze drzewo: ")
        else:
            tree_num = raw_input("Podaj drugie drzewo: ")
        try:
            tree_num = int(tree_num)
        except ValueError:
            if tree_num == 'X':
                break
            else:
                print("Niepoprawny numer drzewa lub znak. Spróbuj ponownie...")
                continue
        if t1 == None:
            t1 = treeFiles[tree_num - 1].eteTree

            print treeFiles[tree_num - 1].eteTree
        else:
            t2 = treeFiles[tree_num - 1].eteTree

            print treeFiles[tree_num - 1].eteTree

    rf = t1.robinson_foulds(t2, unrooted_trees=True)

    print "Dystans RF to "+str(rf[0])+" spośród "+str(rf[0])



class MyTree:
    fileName = ""
    tree = ""
    clusters = ""

    def __contains__(self, root, clade):
        stringElements = clade
        stringElements = stringElements.replace("{", "")
        stringElements = stringElements.split(",")
        for i in range(0, len(stringElements)):
            if stringElements[i] not in root.name:
                return False
        return True

    def __init__(self):
        self.fileName = ""
        self.tree = ""

    def loadFile(self):
        self.fileName = raw_input("Podaj nazwe pliku z drzewem\n")

        try:
            self.tree = Phylo.read(str(self.fileName), "newick")
            print("Wczytano plik " + self.fileName)
        except:
            self.tree = ""
            self.fileName = ""
            print("Wczytywanie nie powiodło się!!!")



def split():
    print("Konwersja drzewa do rodziny zgodnych rozbić\n")

    tree_num = '-1'

    while (True):
        tree_num = raw_input("Podaj drzewo: ")
        try:
            tree_num = int(tree_num)
        except ValueError:
            print("Niepoprawny numer drzewa lub znak. Spróbuj ponownie...")

            continue

        if (tree_num - 1 >= len(treeFiles) or tree_num < 1):
            print("Błąd: Nie ma takiego drzewa.")
        else:
            break

    print("Podane drzewo: ")

    print treeFiles[tree_num - 1].eteTree, "\n\n"

    splitTable = []

    nodes = []

    for node1 in treeFiles[tree_num - 1].eteTree:
        nodes.append(node1.name)
        for node2 in treeFiles[tree_num - 1].eteTree:

            modifiedTree = copy.deepcopy(treeFiles[tree_num - 1].eteTree)

            ancestor = modifiedTree.get_common_ancestor(node1.name, node2.name)
            try:
                modifiedTree.set_outgroup(ancestor)
            except:
                continue


            tableToAdd = []
            for subtree in modifiedTree.children:
                tempTable = []
                for node in subtree:
                    tempTable.append(node.name)
                tableToAdd.append(tempTable)

            if tableToAdd not in splitTable:
                splitTable.append(tableToAdd)

    index = 0

    for node in nodes:
        tableToAdd = [[], []]
        for idx, node2 in enumerate(nodes):
            if index == idx:
                tableToAdd[0].append(node2)
            else:
                tableToAdd[1].append(node2)
        splitTable.append(tableToAdd)
        index += 1

    for table in splitTable:
        print table
    print "\n"

MenuCh = 'x'
readTree()
os.system("clear")
myTree = MyTree()
while MenuCh != '0':
    menu()
    MenuCh = raw_input("Wybierz opcję: ")

    if MenuCh == '1':
        os.system("clear")
        printTrees()

    elif MenuCh == '2':
        os.system("clear")
        consensus()

    elif MenuCh == '3':
        os.system("clear")
        distance()

    elif MenuCh == '4':
        os.system("clear")
        cutLeaf()

    elif MenuCh == '5':
        os.system("clear")
        split()

    elif MenuCh == '0':
        break

    else:
        os.system("clear")
        print("Niepoprawna opcja, spróbuj jeszcze raz")

os.system("clear")
print("Koniec")

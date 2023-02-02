import os
from Bio import Entrez, SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib.pyplot as plt

# Definir la función fasta_downloader
def fasta_downloader():
    """
    Esta función carga id_coati.txt, después se descarga la información correspondiente a los identificadores de acceso
    utilizando ENTREZ de Biopython y, se lo guarda en coati y coati.gb.
    """
    # Cargar id_coati.txt
    with open("id_coati.txt") as file:
        ids = file.read().splitlines()

    # Descargarmos la información utilizando ENTREZ de Biopython
    Entrez.email = "josue.chango@est.ikiam.edu.ec"
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
    records = list(SeqIO.parse(handle, "gb"))
    with open("coati.gb", "w") as output_handle:
        SeqIO.write(records, output_handle, "genbank")
    with open("coati", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

# Definimos la función alignment
def alignment():
    """
    Esta función extrae solo las secuencias de la variable coati y realiza una alineación utilizando ClustalW.
    El resultado se debe guardar en los siguientes archivos los cuales son: coati.aln y coati.dnd en su carpeta de trabajo.
    """
    # Realizar la alineación utilizando ClustalW
    clustalw_cline = ClustalwCommandline("C:\ProgramData\Microsoft\Windows\Start Menu\Programs\ClustalW2\clustalw.exe", infile="coati")
    clustalw_cline()
    os.rename("coati.aln", "coati_aligned.aln")
    os.rename("coati.dnd", "coati_aligned.dnd")

# Definimos la función tree
def tree():
    """
    Esta función ayuda a calcular las distancias utilizando coati.aln y, finalmente imprimimos el árbol filogenético y
    lo guardamos como coati_phylotree.pdf en nuestro repositorio
    """
    # Cargar las secuencias que ya están alineadas
    with open("coati_aligned.aln", "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    
    # Calculamos las distancias
    calculator = DistanceCalculator("blosum62")
    dm = calculator.get_distance(records)

    # Construir el árbol
    constructor = Distance

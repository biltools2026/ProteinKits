ProteinKits is a program for protein structure prediction.

Some running program are as below:

Example 1
perl run.pl -3d E:\Udisk128\ProteinKits\fold\example\seq.fasta

The example predict the 3D model of the sequence using local database

Example 2
perl run.pl -3d E:\Udisk128\ProteinKits\fold\example\seq.fasta -webpdb E:\Udisk128\ProteinKits\fold\example\seq.fasta

The example predict the 3D model of the sequence using local database and online PDB database

Example 3
perl run.pl -3d E:\Udisk128\ProteinKits\fold\8tim\seq.fasta -t E:\Udisk128\ProteinKits\fold\PDB\8tim.pdb -webpdb E:\Udisk128\ProteinKits\fold\8tim\seq.fasta

The example predict the 3D model of the sequence using the template 8tim.pdb, local database and online PDB database






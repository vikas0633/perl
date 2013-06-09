''-----------------------------------------------------------+
'' batch_ltrstruc.vbs - Run LTR_Struc in batch mode          |
''-----------------------------------------------------------+
''                                                           |
''  AUTHOR: James C. Estill                                  |
'' CONTACT: JamesEstill_at_gmail.com                         |
'' STARTED: 12/13/2006                                       |
'' DESCRIPTION:                                              |
''  Faciliates running the LTR_Struc program in              |
''  a batch mode. This prevents the user from needing        |
''  to answer command line prompts and hit the               |
''  enter key for each BAC. By default this program will     |
''  run LTR_Struc in level 1 to find the most LTR Retros.    |
''                                                           |
'' REQUIRED SOFTWARE:                                        |
''  LTR_Struc can be obtained at:                            |
''  http://www.genetics.uga.edu/retrolab/data/LTR_Struc.html |
''                                                           |
'' SEE ALSO:                                                 |
''  DAWG-PAWS                                                |
''  This program is part of the DAWG-PAWS  set of            |
''  scripts designed for the annotation of plant             |
''  genomes. DAWG-PAWS is available from SourceForge         |
''  availalbe from SourceForge.                              |
''  http://dawgpaws.sourceforge.net/                         |
''                                                           |
''                                                           |
'' From the command line this can be launched                |
'' with the Cscript program                                  |
'' Cscript.exe //NoLogo batch_ltrstruc.vbs | cmd.exe         |
''-----------------------------------------------------------+

Wscript.Stdout.WriteLine "LTR_STRUC_1_1"
Wscript.Sleep 1000
Wscript.Stdout.WriteLine "Y"
Wscript.Sleep 1000
'' This select N for using the standard settings
Wscript.Stdout.WriteLine "N"
'' This selects 1 for the most stringent
Wscript.Stdout.WriteLine "1"
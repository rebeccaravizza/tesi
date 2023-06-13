########################################
	INFORMAZIONI UTILIZZO
########################################

All'interno del codice identification.m ci sono 2 importanti variabili da verificare:

- parentFolder: Deve contenere il path alla cartella dove sono contenute tutte le cartelle con i dati preprocessati (ex. preprocessedData)
- savePath: Deve contenere il path dove vogliono essere salvati i dati (può essere lasciato anche quello di default)

########################################
	UTILIZZO execution_registry
########################################

Automaticamente il programma scriverà nel file execution_registry i nomi dei file che sono stati processati,
se si vogliono escludere dei file da eseguire, basterà scriverli lì, allo stesso modo, se si vuole riprocessare un file
sarà necessario cancellarlo manualmente da li.
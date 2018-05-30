## Aplicacion para graficar las funciones de transferencia de todo el dominio 
y los diferencias entre secciones rectas y curvas: 


Al inicio del codigo se le deben definir la ruta raiz del caso canonico es decir,
la ruta donde estaran los datos y las graficas. En este trabajo en las carpetas:
01_Canonicos a 06_Canonicos. En cada caso se deben definir tres forder internos

Folder: 1_Dates

Contiene los datos de las funciones de transferencia para los movimientos P,SV,SH en fortamo .txt. 
En total se deben especificar 10 archivos 

FTX_A_CSV, FTX_A_RSV, FTX_A_CSV, FTY_A_RSV, FTY_A_CSV, FTY_A_RSV, FTY_A_CSV, FTY_A_RSV,FTZ_RSH, FTZ_CSH. 

Tambien es preciso definir la variable Type para el tipo de canon (C o T) y la variable canon.  Ademas definir el archivo dates.py

collection: contiene los Ld analizados estos es para hacer todos los valores de Ld para una geometría seleccionada.
En caso que se corran distintos valores distintos valores de Ld se puede activar la opcion de recorrer las carpetas
con la función fn.Setpath(Path,NL) en linea 90, sino se apaga  y se debe dejar prendida la asignacion para Path1 en
linea 91.


Los resultados se entregan: 

Folder: 2_Error
Graficas del error asociado a cada frecuencia y por cada tipo de onda. 

ej: Error_P_0.5. Es el error para onda P con frecuencia de 0.5 hz

Folder: 3_FT
Graficas de funcion de transferencia para ondas P,SV,SH para una triada frecuencia

ej: FT_Tot_0: Función de transferencia para la primera triada 



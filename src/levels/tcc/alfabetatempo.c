/* set yrange [1.0:3.0] */
/* set ylabel "beta" */
/* set xrange [0.5:4.0] */
/* set xlabel "alfa" */
/* set pm3d */
/* splot "param.dat" with lines  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libpq-fe.h>
#include <math.h>

void
do_exit (PGconn *conn) {
  PQfinish(conn);
  exit(1);
}

int
main () {
  int i;
  int nparamdb = 2;
  const char *paramValues[2];
  
  /* FILE *arq = fopen ("param.dat", "w"); */

  for (i = 0; i < nparamdb; ++i) 
    paramValues[i] = (char*) malloc (4*sizeof (char));

  PGconn *conn = PQconnectdb("user=allisson dbname=acodb");

  if (PQstatus(conn) == CONNECTION_BAD){
    printf ("connection to database failsed\n");
    do_exit (conn);
  }
  double alfa, beta;
  int somalen = 0;
  float medialen = 0;
  float desviolen = 0;
  int maiorlen = 0;
  int desv1 = 0;
  int exec, execucoes = 30;

  PGresult *res;

  //BEGIN SELECT MAX LENGTH
  /* char *s = "SELECT MAX(len) FROM antsolutions"; */
  /* res = PQexec (conn, s); */

  /* if (PQresultStatus (res) != 2) { */
  /*   printf ("No data retrieved\n"); */
  /*   PQclear(res); */
  /*   do_exit(conn); */
  /* } */
  /* /\* maiorlen = 172; *\/ */
  /* maiorlen = atoi (PQgetvalue (res, 0, 0)); */
  /* printf ("maiorlen = %d\n",maiorlen); */
  /* PQclear (res); */
  //END SELECT MAX LENGTH

  for (alfa = 0.5; alfa <= 3.5; alfa+=0.5) {
    for (beta = 0.5; beta <= 4.; beta+=0.5) {

      /* for (exec = 0; exec < execucoes; ++exec) { */
	sprintf (paramValues[0], "%.1lf", alfa);
	sprintf (paramValues[1], "%.1lf", beta);
	/* sprintf (paramValues[2], "%d", exec);   */
	char *stm="SELECT AVG(timesol) FROM antsolutions WHERE alfa=$1 AND beta=$2";
	res=PQexecParams(conn, stm, 2, NULL, paramValues, NULL, NULL, 0);
      
	if (PQresultStatus (res) != PGRES_TUPLES_OK) {
	  printf ("No data retrieved\n");
	  PQclear(res);
	  do_exit(conn);
	}

	int rows = PQntuples (res);
	int leng = 0;
	double tempo;
	for (i = 0; i < rows; ++i) {

	  tempo = atof (PQgetvalue (res, i ,0));
	  /* somalen += atoi (PQgetvalue (res, i ,10)); */
	  /* printf ("alfa = %s\tbeta = %s\n", */
	  /* 	    PQgetvalue (res, i, 4),PQgetvalue (res, i, 5)); */
	}
	/* printf ("exec %d -> len = %d\n", exec, leng); */
	/* printf ("rows = %d\t somalen = %d\n", rows, somalen); */
	/* printf ("Para alfa = %s e beta = %s\na média de len é %d\n\n", PQgetvalue (res, 0, 4), PQgetvalue (res, 0, 5), somalen/rows); */
	PQclear (res);
      /* } //markexec */

            printf ("alfa = %.1lf, beta = %.1lf, tempo = %3.2f\n", alfa, beta, tempo);
      /* printf ("alfa = %.1lf, beta = %.1lf, len = %3.2f, desvio = %3.2f qualidade = %3.2lf\n", alfa, beta, medialen, desviolen,  maiorlen - medialen); */

      /* fprintf (arq, "%.1lf %.1lf %.1lf\n", alfa, beta, maiorlen - medialen); */
    } //markbeta
    putchar ('\n');
    /* fprintf (arq, "\n"); */
  } //markalfa

  PQfinish (conn);
  /* fclose (arq); */
  
  return 0;
}
  

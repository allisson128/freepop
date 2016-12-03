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
  int nparamdb = 4;
  const char *paramValues[4];
  
  FILE *arq100 = fopen ("new100formigas.dat", "w");
  /* FILE *arq90 = fopen ("new90formigas.dat", "w"); */
  /* FILE *arq80 = fopen ("new80formigas.dat", "w"); */
  /* FILE *arq70 = fopen ("new70formigas.dat", "w"); */
  
  for (i = 0; i < nparamdb; ++i) 
    paramValues[i] = (char*) malloc (4*sizeof (char));

  PGconn *conn = PQconnectdb("user=allisson dbname=acodb");

  if (PQstatus(conn) == CONNECTION_BAD){
    printf ("connection to database failsed\n");
    do_exit (conn);
  }
  double alfa, beta;
  int somalen = 0;
  int medialen = 0;
  float desviolen = 0;
  int maiorlen = 0;
  int desv1 = 0;
  int exec, execucoes = 30;
  int step, steps = 10000;
  int ants;
  
  PGresult *res;
  /* int s, o, n, c; */
  /* s = o = n = c = 1; */
  
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

  for (i = 0; i < 300; ++i)
    fprintf (arq100, "%d %d\n", i, 0);

  for (alfa = 2.5; alfa <= 2.51; alfa+=0.5) {
    for (beta = 1.5; beta <= 1.51; beta+=0.5) {

      for (ants = 100; ants <= 100; ants+= 10) {

	/* for (step = 0; step < steps; ++step) { */
	step = 300;
	while (step < steps) {

	sprintf (paramValues[0], "%.1lf", alfa);
	sprintf (paramValues[1], "%.1lf", beta);
	sprintf (paramValues[2], "%d", ants);
	sprintf (paramValues[3], "%d", step);
	char *stm="SELECT COUNT(*) FROM antsolutions WHERE len <= 71 AND alfa=$1 AND beta=$2 AND ants = $3 AND step <= $4";
	res=PQexecParams(conn, stm, 4, NULL, paramValues, NULL, NULL, 0);
      
	if (PQresultStatus (res) != PGRES_TUPLES_OK) {
	  printf ("No data retrieved\n");
	  PQclear(res);
	  do_exit(conn);
	}

	int rows = PQntuples (res);
	int leng = 0;
	for (i = 0; i < rows; ++i) {
	  medialen  = atoi (PQgetvalue (res, i ,0)) / execucoes;
	}


	printf("ants = %d, step = %5d, freq = %5d\n",ants,step, medialen);

	/* if (ants == 70) */
	/*   fprintf (arq70, "%d %d\n", step, medialen); */
	/* else if (ants == 80) { */
	/*   fprintf (arq80, "%d %d\n", step, medialen); */
	/*   if (s) { */
	/*     fclose (arq70); */
	/*     s = 0; */
	/*   } */
	/* } */
	/* else if (ants == 90) { */
	/*   fprintf (arq90, "%d %d\n", step, medialen); */
	/*   if (o) { */
	/*     fclose (arq80); */
	/*     o = 0; */
	/*   } */
	/* } */
	/* else if (ants == 100) { */
	  fprintf (arq100, "%d %d\n", step, medialen);
	  /* if (n) { */
	  /*   fclose (arq90); */
	  /*   n = 0; */
	  /* } */
	/* } */
	PQclear (res);
	  
	if (step < 1000)
	  step += 1;
	else if (step < 5000)
	  step += 10;
	else
	  step += 100;
      } // markstep
      } //markant

    } //markbeta
    /* putchar ('\n'); */
    /* fprintf (arq, "\n"); */
  } //markalfa

  PQfinish (conn);


  fclose (arq100);
  
  return 0;
}
  

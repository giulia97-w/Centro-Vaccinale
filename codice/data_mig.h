 

//tempi sistema
#define STOP   300.0               
#define INT_MAX      2147483647
//numero di servers
#define SERVERS_ACCETTAZIONE 2
#define SERVERS_ANAMNESI 2
#define SERVERS_INOCULAZIONE 2
#define SERVERS_INOCULAZIONE_ATTENUATO 1
//lag autocorrelazione
#define LAG 2
//tasso di arrivo
#define LAMBDA    0.456
//probabilità
#define prob_no_vaccino 1.0//prob di chi può fare il vaccino quindi non esonerato da esso
#define prob_vaccino_attenuato 0.23
//Batch e K
#define B 4096
#define K 64
#define SEED 123456789
#define ITERATIONS 64
#define LOC              0.95 
//tassi di servizio
#define MU_ACCETTAZIONE 0.515
#define MU_ANAMNESI  0.3
#define MU_INOCULAZIONE  0.249
#define MU_INOCULAZIONE_ATTENUATO 0.36
//numero eventi
#define EVENTI 9
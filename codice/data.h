 

//tempi sistema
#define STOP  360.0               
#define INT_MAX      2147483647
//numero di servers
#define SERVERS_ACCETTAZIONE 1
#define SERVERS_ANAMNESI 3
#define SERVERS_INOCULAZIONE 2
#define SERVERS_INOCULAZIONE_ATTENUATO 2
//lag autocorrelazione
#define LAG 1
//tasso di arrivo
#define LAMBDA    0.456
//probabilità
#define prob_no_vaccino 1.0//prob di chi può fare il vaccino quindi non esonerato da esso
#define prob_vaccino_attenuato 0.218978
//Batch e K
#define B 2048
#define K 128
#define SEED 123456789
#define ITERATIONS 128
#define LOC              0.95 
//tassi di servizio
#define MU_ACCETTAZIONE 0.515
#define MU_ANAMNESI  0.3
#define MU_INOCULAZIONE  0.249
#define MU_INOCULAZIONE_ATTENUATO 0.36
//numero eventi
#define EVENTI 9
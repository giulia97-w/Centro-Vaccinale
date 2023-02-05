typedef struct multiserver{

    int busy;
    double service;
    long served;


}multiserver;

typedef struct {
    double current;
    double x;
}event_list;

double sim();
void PlantSeeds(long x);

typedef struct center{
    char *name;
    double area;
    double service;
    double index;
    double servers;
    double init_time;
    double lastArrival;
    double current;
    double queue;
    double number;
    double mean;
    
    
    

}center;



struct t{
   double current;
   double next;
   double last; 

}t;

typedef struct dati{
    double lambda;
    double Ts;
    double Tq;
    double Es;
    double N;
    double Nq;
    double ro;
    char *name;
    double job;

}dati;



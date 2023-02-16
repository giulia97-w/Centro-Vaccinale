#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rvms.h"
#include "config.h"
#include "rngs.h"
#include "rvgs.h"
#include "data.h"



double arrival = 0.0;
double currentBatch = 0.0;


unsigned long long fattoriale(unsigned int n) {
    unsigned long long result = 1;
    for (unsigned int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}


int findServer(multiserver multiserver[], int off_num){
    int i, num, lowest;

    num = -1;
    lowest = INT_MAX;

    for(i=0; i<off_num; i++){
        if(multiserver[i].busy == 0 && multiserver[i].service < lowest){
            lowest = multiserver[i].service;
            num = i;
        }
    }

    return num;
}

int batch(int a, int n){

    if(a<n) {
      return 1;  
    }else{
        return 0;
    }

}

double GetArrival(){
    
  SelectStream(6);
  arrival += Exponential(1/LAMBDA);
  return (arrival);
}


double GetService( char *name){
    if(strcmp(name, "accettazione") == 0){
        SelectStream(1);
        return Exponential(1/MU_ACCETTAZIONE);
    }
    else if(strcmp(name, "anamnesi")==0){
        SelectStream(2);
        return Exponential(1/MU_ANAMNESI);
    }else if(strcmp(name, "inoculazione")==0){
        SelectStream(3);
        return Exponential(1/MU_INOCULAZIONE);

    }else {
        SelectStream(4);
        return Exponential(1/MU_INOCULAZIONE_ATTENUATO);
    }
}


int NextEvent(event_list events[]){
  int e;                                      
  int i = 0;

  while (events[i].x == 0)       /* find the index of the first 'active' */
    i++;                        /* element in the event list            */ 
  e = i;                        
  while (i < EVENTI) {         /* now, check the others to find which  */
    i++;                        /* event type is most imminent          */
    if ((events[i].x == 1) && (events[i].current < events[e].current))
      e = i;
  }
  return (e);
}



int batches[5]; 
double STOP1;
double n;
double sim(char *simulation, dati datiFinite[], dati matrix[B][4]){
    
    

    multiserver sumAnamnesi[SERVERS_ANAMNESI];
    multiserver sumInoculazione[SERVERS_INOCULAZIONE];
    multiserver sumInoculazione_attenuato[SERVERS_INOCULAZIONE_ATTENUATO];
    
//Reset arrays for servers
int servers[] = {SERVERS_ANAMNESI, SERVERS_INOCULAZIONE, SERVERS_INOCULAZIONE_ATTENUATO};
for(int s=0; s<3; s++){
    for(int i=0; i<servers[s]; i++){
        if(s==0) sumAnamnesi[i].served = 0, sumAnamnesi[i].busy = 0, sumAnamnesi[i].service = 0.0;
        if(s==1) sumInoculazione[i].served = 0, sumInoculazione[i].busy = 0, sumInoculazione[i].service = 0.0;
        if(s==2) sumInoculazione_attenuato[i].served = 0, sumInoculazione_attenuato[i].busy = 0, sumInoculazione_attenuato[i].service = 0.0;
    }
}

    event_list fromAnamnesi[SERVERS_ANAMNESI];
    event_list fromInoculazione[SERVERS_INOCULAZIONE];
    event_list fromInoculazione_attenuato[SERVERS_INOCULAZIONE_ATTENUATO];
    event_list event[EVENTI];


int numCenters = 5;
    for(int i = 0; i< numCenters; i++){
        batches[i] = 0;
    }

    for(int s = 0; s < 3; s++) {
        for (int i = 0; i < servers[s]; i++) {
            if(s==0) fromAnamnesi[i].x = 0, fromAnamnesi[i].current = INT_MAX;
            if(s==1) fromInoculazione[i].x = 0, fromInoculazione[i].current = INT_MAX;
            if(s==2) fromInoculazione_attenuato[i].x = 0, fromInoculazione_attenuato[i].current = INT_MAX;
        }
    }

    for (int i = 0; i < EVENTI; i++) {
        event[i].x = 0;
        event[i].current = INT_MAX;
    }

    center accettazione;
    center anamnesi;
    center inoculazione;
    center inoculazione_attenuato;
    center centro_vaccinale;

    accettazione.queue = 0.0;
    accettazione.service = 0.0;
    accettazione.index = 0.0;
    accettazione.area = 0.0;
    accettazione.name = "accettazione";
    accettazione.init_time= accettazione.lastArrival;
    accettazione.lastArrival = 0.0; 
    accettazione.current = 0.0;
    accettazione.servers = SERVERS_ACCETTAZIONE;
    accettazione.number = 0.0;

    anamnesi.queue = 0.0;
    anamnesi.service = 0.0;
    anamnesi.index = 0.0;
    anamnesi.area = 0.0;
    anamnesi.name = "anamnesi";
    anamnesi.init_time= anamnesi.lastArrival;
    anamnesi.lastArrival = 0.0;
    anamnesi.current = 0.0;
    anamnesi.servers = SERVERS_ANAMNESI;
    anamnesi.number = 0.0;

    inoculazione.queue = 0.0;
    inoculazione.service = 0.0;
    inoculazione.index = 0.0;
    inoculazione.area = 0.0;
    inoculazione.name = "inoculazione";
    inoculazione.init_time= inoculazione.lastArrival;
    inoculazione.lastArrival = 0.0;
    inoculazione.current = 0.0;
    inoculazione.servers = SERVERS_INOCULAZIONE;
    inoculazione.number = 0.0;

    inoculazione_attenuato.queue = 0.0;
    inoculazione_attenuato.service = 0.0;
    inoculazione_attenuato.index = 0.0;
    inoculazione_attenuato.area = 0.0;
    inoculazione_attenuato.name = "inoculazione_attenuato";
    inoculazione_attenuato.init_time= inoculazione_attenuato.lastArrival;
    inoculazione_attenuato.lastArrival = 0.0;
    inoculazione_attenuato.current = 0.0;
    inoculazione_attenuato.servers = SERVERS_INOCULAZIONE_ATTENUATO;
    inoculazione_attenuato.number = 0.0;



    arrival = 0.0;
    
    
    

    event[0].current = GetArrival();
    event[0].x = 1;
    double inBatch = event[0].current;

    

    while(((strcmp(simulation,"finite")==0) && ((event[0].current < STOP) || ((accettazione.number + anamnesi.number + inoculazione.number + inoculazione_attenuato.number) > 0))) || ((strcmp(simulation,"infinite")==0) && (currentBatch < K * 4))){
        
        int e = NextEvent(event);
        t.next = event[e].current;


//integrali accettazione
    if(accettazione.number > 0){
        accettazione.area += (t.next - t.current) * accettazione.number;
        
            if(accettazione.number > accettazione.servers) {
                
                accettazione.service += (t.next - t.current)*accettazione.servers;
            }
            else{
                accettazione.service +=(t.next - t.current)*accettazione.number;
            }
            
            if(accettazione.number > accettazione.servers - 1){
                accettazione.queue += (t.next - t.current) * (accettazione.number - accettazione.servers);

                }    }


//integrali anamnesi
    if(anamnesi.number > 0){
        anamnesi.area += (t.next - t.current) * anamnesi.number;
        
            if(anamnesi.number > anamnesi.servers) {
                
                anamnesi.service += (t.next - t.current)*anamnesi.servers;
            }
            else{
                anamnesi.service +=(t.next - t.current)*anamnesi.number;
            }
            
            if(anamnesi.number > anamnesi.servers - 1){
                anamnesi.queue += (t.next - t.current) * (anamnesi.number - anamnesi.servers);

                }    }


//integrali inoculazione        
    if(inoculazione.number>0){
        inoculazione.area += (t.next - t.current) * inoculazione.number;
       
            if(inoculazione.number > inoculazione.servers) {
                
                inoculazione.service += (t.next - t.current)*inoculazione.servers;
            }
            else{
                inoculazione.service +=(t.next - t.current)*inoculazione.number;
            }
            
            if(inoculazione.number > inoculazione.servers){
                inoculazione.queue += (t.next - t.current) * (inoculazione.number - inoculazione.servers);

            } }


//integrali inoculazione attenuato
    if(inoculazione_attenuato.number>0){
        inoculazione_attenuato.area += (t.next - t.current) * inoculazione_attenuato.number;
       
            if(inoculazione_attenuato.number > inoculazione_attenuato.servers) {
                
                inoculazione_attenuato.service += (t.next - t.current)*inoculazione_attenuato.servers;
            }
            else{
                inoculazione_attenuato.service +=(t.next - t.current)*inoculazione_attenuato.number;
            }
            
            if(inoculazione_attenuato.number > inoculazione_attenuato.servers - 1){
                inoculazione_attenuato.queue += (t.next - t.current) * (inoculazione_attenuato.number - inoculazione_attenuato.servers);

            }
    }


        t.current = t.next;    

            
        if(e == 0){
        

                centro_vaccinale.number++;
                accettazione.number++; 
                accettazione.lastArrival = t.current;
  
                event[0].current = GetArrival();

                if(strcmp(simulation, "finite")== 0){
                if(event[0].current> STOP){
                    t.last = t.current;
                    event[0].x = 0;
                    event[0].current = INT_MAX;
                }
            }
            
                if(accettazione.number == 1){
                    
                    event[1].current = t.current + GetService("accettazione");
                    event[1].x = 1;
                }
            
            }
        

        


        if(e == 1){

            accettazione.index++;
            accettazione.number--;
            accettazione.current = t.current;
            
            
            if(accettazione.number > 0){
                
                event[1].current = t.current + GetService("accettazione");
                event[1].x = 1;
            }else{
                event[1].current = INT_MAX;
            }
            
                event[2].current = event[1].current;
                event[2].x = 1;
                
                e = 2; 
                
        }
        if(e == 2){
            
            double fromcenter = INT_MAX;
            anamnesi.number++;
            anamnesi.lastArrival = t.current;
            
            event[2].x = 0;

            
            if(anamnesi.number <= SERVERS_ANAMNESI){
                
                double service = GetService("anamnesi");
                int s = findServer(sumAnamnesi, SERVERS_ANAMNESI);

                sumAnamnesi[s].service += service;
                sumAnamnesi[s].served++;
                sumAnamnesi[s].busy = 1;
                fromAnamnesi[s].current = t.current + service;
                fromAnamnesi[s].x = 1;
                

                event[3].x = 1;
                for (int i = 0; i<SERVERS_ANAMNESI; i++){
                        if(fromAnamnesi[i].x == 1){
                            if (fromAnamnesi[i].current < fromcenter){
                                fromcenter = fromAnamnesi[i].current;
                            }
                        } 
                            
                        }
                
                
                event[3].current = fromcenter;
                
        
            }
            }

        if(e == 3){  
            
           
            double p = Uniform(0,1);
            double p1 = Uniform(0,1);
            double fromcenter = INT_MAX;
            anamnesi.index++;
            anamnesi.number--;
            anamnesi.current = t.current;
            
            

            int serverIndex1 = 0;

            for(int i = 0; i<SERVERS_ANAMNESI; i++){
                if(fromAnamnesi[i].current == event[3].current){
                    serverIndex1 = i;
                }
            }
            
            if(anamnesi.number >=SERVERS_ANAMNESI){
                
                double service = GetService("anamnesi");

                sumAnamnesi[serverIndex1].service += service;
                sumAnamnesi[serverIndex1].served++;
                fromAnamnesi[serverIndex1].current = t.current + service;
                    
            }else{
                sumAnamnesi[serverIndex1].busy = 0;
                fromAnamnesi[serverIndex1].x = 0;
                fromAnamnesi[serverIndex1].current = INT_MAX;

            }
            for (int i = 0; i<SERVERS_ANAMNESI; i++){
                        if(fromAnamnesi[i].x == 1){
                            if(fromAnamnesi[i].current < fromcenter){
                                fromcenter = fromAnamnesi[i].current;
                            
                            }
                        } 
                            
                        }
              
            event[3].current =  fromcenter;
          
                
            
            SelectStream(7);
        if(p1 < prob_no_vaccino){
            if(p<(prob_vaccino_attenuato)){
                
                
                event[4].current = event[3].current;
                    event[4].x = 1;
                
                    e = 4;
            }else{
                SelectStream(8);
                
                event[6].current = event[3].current;
                event[6].x = 1;
                        
                e = 6;
            }
                  
            }}

        if (e == 4){
            
           
          
            double fromcenter = INT_MAX;
            
            
            inoculazione_attenuato.number++;
            inoculazione_attenuato.lastArrival = t.current;
            
            event[4].x = 0;

            
            if(inoculazione_attenuato.number <= SERVERS_INOCULAZIONE_ATTENUATO){
                
                double service = GetService("inoculazione_attenuato");

                int s = findServer(sumInoculazione_attenuato, SERVERS_INOCULAZIONE_ATTENUATO);

                sumInoculazione_attenuato[s].service += service;
                sumInoculazione_attenuato[s].served++;
                sumInoculazione_attenuato[s].busy = 1;
                fromInoculazione_attenuato[s].current = t.current + service;
                fromInoculazione_attenuato[s].x = 1;
                

                event[5].x = 1;
                for (int i = 0; i<SERVERS_INOCULAZIONE_ATTENUATO; i++){
                        if(fromInoculazione_attenuato[i].x == 1){
                            if (fromInoculazione_attenuato[i].current < fromcenter){
                                fromcenter = fromInoculazione_attenuato[i].current;
                            }
                        } 
                            
                        }
                    
                event[5].current = fromcenter;
                
        

        }}

        if(e == 5){
            double fromcenter = INT_MAX;
            inoculazione_attenuato.index++;
            inoculazione_attenuato.number--;
            inoculazione_attenuato.current = t.current;
            
            centro_vaccinale.index++;
            centro_vaccinale.number--;
            

            int serverIndex1 = 0;

            for(int k = 0; k<SERVERS_INOCULAZIONE_ATTENUATO; k++){
                if(fromInoculazione_attenuato[k].current == event[5].current){
                    serverIndex1 = k;
                }
            }
            
            if(inoculazione_attenuato.number >=SERVERS_INOCULAZIONE_ATTENUATO){
                
                double service = GetService("inoculazione_attenuato");

                sumInoculazione_attenuato[serverIndex1].service += service;
                sumInoculazione_attenuato[serverIndex1].served++;
                fromInoculazione_attenuato[serverIndex1].current = t.current + service;
                    
            }else{
                sumInoculazione_attenuato[serverIndex1].busy = 0;
                fromInoculazione_attenuato[serverIndex1].x = 0;
                fromInoculazione_attenuato[serverIndex1].current = INT_MAX;

            }
            for (int i = 0; i<SERVERS_INOCULAZIONE_ATTENUATO; i++){
                        if(fromInoculazione_attenuato[i].x == 1){
                            if(fromInoculazione_attenuato[i].current < fromcenter){
                                fromcenter = fromInoculazione_attenuato[i].current;
                            }
                        } 
                            
                        }
                    
            event[5].current =  fromcenter;
            

        }
        
        if(e == 6){
            
            double fromcenter = INT_MAX;
            inoculazione.number++;
            inoculazione.lastArrival = t.current;
            event[6].x = 0;

           
            if(inoculazione.number <= SERVERS_INOCULAZIONE){
                
                double service = GetService("inoculazione");

                int s = findServer(sumInoculazione, SERVERS_INOCULAZIONE);

                sumInoculazione[s].service += service;
                sumInoculazione[s].served++;
                sumInoculazione[s].busy = 1;
                fromInoculazione[s].current = t.current + service;
                fromInoculazione[s].x = 1;
                

                event[7].x = 1;
                for (int i = 0; i<SERVERS_INOCULAZIONE; i++){
                        if(fromInoculazione[i].x == 1){
                            if (fromInoculazione[i].current < fromcenter){
                                fromcenter = fromInoculazione[i].current;
                            }
                        } 
                            
                        }
                    
                event[7].current = fromcenter;
                }
        
        } 

        if(e == 7){    
            
            double fromcenter = INT_MAX;
            inoculazione.index++;
            inoculazione.number--;
           
            inoculazione.current = t.current;
            centro_vaccinale.index++;
            centro_vaccinale.number--;
            

            int serverIndex1 = 0;

            for(int k = 0; k<SERVERS_INOCULAZIONE; k++){
                if(fromInoculazione[k].current == event[7].current){
                    serverIndex1 = k;
                }
            }
            
            if(inoculazione.number >=SERVERS_INOCULAZIONE){
                
                double service = GetService("inoculazione");

                sumInoculazione[serverIndex1].service += service;
                sumInoculazione[serverIndex1].served++;
                fromInoculazione[serverIndex1].current = t.current + service;
                    
            }else{
                sumInoculazione[serverIndex1].busy = 0;
                fromInoculazione[serverIndex1].x = 0;
                fromInoculazione[serverIndex1].current = INT_MAX;

            }
            for (int i = 0; i<SERVERS_INOCULAZIONE; i++){
                        if(fromInoculazione[i].x == 1){
                            if(fromInoculazione[i].current < fromcenter){
                                fromcenter = fromInoculazione[i].current;
                            }
                        } 
                            
                        }
                    
            event[7].current =  fromcenter;
            
                
        } 
        

    if(strcmp(simulation, "infinite") == 0){  
        
            
                if((accettazione.index == B) && (batches[0] < K)){
                
                        dati out;
                        out.name = accettazione.name;               
                        out.lambda = accettazione.index / (accettazione.lastArrival - accettazione.init_time);
                        out.job = accettazione.index;
                        out.Ts = accettazione.area / accettazione.index;
                        out.Tq = accettazione.queue / accettazione.index;
                        out.Es = accettazione.service/ accettazione.index;
                        out.N = (accettazione.area) / (accettazione.current - accettazione.init_time);
                        out.Nq = accettazione.queue / (accettazione.current - accettazione.init_time);
                        out.ro = accettazione.service /(accettazione.current - accettazione.init_time);
                        
                        
                        currentBatch++;
                        matrix[batches[0]][0] = out;
                        batches[0] = batches[0] + 1;
                        
                        accettazione.queue = 0.0;
                        accettazione.service = 0.0;
                        accettazione.index = 0.0;
                        accettazione.area = 0.0;
                        accettazione.name = "accettazione";
                        accettazione.init_time= accettazione.lastArrival;
                        accettazione.lastArrival = 0.0;
                        accettazione.current = 0.0;
                        accettazione.servers = SERVERS_ACCETTAZIONE;
             }           
            
        
            if((anamnesi.index == B) && (batches[1] < K)){

                dati out;
                out.name = anamnesi.name;
                out.lambda = (anamnesi.index / (anamnesi.lastArrival - anamnesi.init_time));
                out.job = anamnesi.index;
                out.Ts = (anamnesi.area) / anamnesi.index;
                out.Tq = anamnesi.queue / anamnesi.index;
                out.Es = (anamnesi.service/anamnesi.servers) / anamnesi.index;
                out.N = anamnesi.area / (anamnesi.current - anamnesi.init_time);
                out.Nq = anamnesi.queue / (anamnesi.current - anamnesi.init_time);
                out.ro = anamnesi.service / (anamnesi.servers*(anamnesi.current - anamnesi.init_time));
                            
                currentBatch++;
                matrix[batches[1]][1] = out;
                batches[1] = batches[1] +1;

                anamnesi.queue = 0.0;
                anamnesi.service = 0.0;
                anamnesi.index = 0.0;
                anamnesi.area = 0.0;
                anamnesi.name = "anamnesi";
                anamnesi.init_time= anamnesi.lastArrival;
                anamnesi.lastArrival = 0.0;
                anamnesi.current = 0.0;
                anamnesi.servers = SERVERS_ANAMNESI;
                
                
            }
        
            if((inoculazione.index == B) && (batches[2] < K)){
            
                dati out;
                out.name = inoculazione.name;
                out.lambda = inoculazione.index / (inoculazione.lastArrival - inoculazione.init_time);
                out.job = inoculazione.index;
                out.Ts = inoculazione.area / inoculazione.index;
                out.Tq = inoculazione.queue / inoculazione.index;
                out.Es = (inoculazione.service/inoculazione.servers) / inoculazione.index;
                out.N = (inoculazione.area) / (inoculazione.current - inoculazione.init_time);
                out.Nq = inoculazione.queue / (inoculazione.current - inoculazione.init_time);
                out.ro = inoculazione.service / (inoculazione.servers*(inoculazione.current - inoculazione.init_time));
                currentBatch++;
                matrix[batches[2]][2] = out; 
                batches[2] = batches[2] + 1;


                inoculazione.queue = 0.0;
                inoculazione.service = 0.0;
                inoculazione.index = 0.0;
                inoculazione.area = 0.0;
                inoculazione.name = "inoculazione";
                inoculazione.init_time= inoculazione.lastArrival;
                inoculazione.lastArrival = 0.0;
                inoculazione.current = 0.0;
                inoculazione.servers = SERVERS_INOCULAZIONE;
                
                
            }
        
            if((inoculazione_attenuato.index == B) && (batches[3] < K)){
            
                dati out;
                out.name = inoculazione_attenuato.name;
                out.lambda = inoculazione_attenuato.index / (inoculazione_attenuato.lastArrival - inoculazione_attenuato.init_time);
                out.job = inoculazione_attenuato.index;
                out.Ts = inoculazione_attenuato.area / inoculazione_attenuato.index;
                out.Tq = inoculazione_attenuato.queue / inoculazione_attenuato.index;
                out.Es = (inoculazione_attenuato.service/inoculazione_attenuato.servers) / inoculazione_attenuato.index;
                out.N = (inoculazione_attenuato.area) / (inoculazione_attenuato.current - inoculazione_attenuato.init_time);
                out.Nq = inoculazione_attenuato.queue / (inoculazione_attenuato.current - inoculazione_attenuato.init_time);
                out.ro = inoculazione_attenuato.service / (inoculazione_attenuato.servers*(inoculazione_attenuato.current - inoculazione_attenuato.init_time));
                currentBatch++;
                matrix[batches[3]][3] = out; 
                batches[3] = batches[3] + 1;


                inoculazione_attenuato.queue = 0.0;
                inoculazione_attenuato.service = 0.0;
                inoculazione_attenuato.index = 0.0;
                inoculazione_attenuato.area = 0.0;
                inoculazione_attenuato.name = "inoculazione_attenuato";
                inoculazione_attenuato.init_time= inoculazione_attenuato.lastArrival;
                inoculazione_attenuato.lastArrival = 0.0;
                inoculazione_attenuato.current = 0.0;
                inoculazione_attenuato.servers = SERVERS_INOCULAZIONE_ATTENUATO;
                
                
            }

           
            if((centro_vaccinale.index == B) && (batches[4] < K)){
            
                double obsTime = centro_vaccinale.lastArrival - inBatch;
                inBatch = t.current;
            

                int batch = batches[4];
           
                batches[4]++;
                centro_vaccinale.queue = 0.0;
                centro_vaccinale.service = 0.0;
                centro_vaccinale.index = 0.0;
                centro_vaccinale.area = 0.0;
                centro_vaccinale.name = "centro_vaccinale";
                centro_vaccinale.init_time= centro_vaccinale.lastArrival;
                centro_vaccinale.lastArrival = 0.0;
                centro_vaccinale.current = 0.0;
                centro_vaccinale.servers = 0;
                centro_vaccinale.number = 0.0;

            }}}
        
    
    


    if((strcmp(simulation,"finite")==0)){

        dati accettazioneOutput;
        accettazioneOutput.name = accettazione.name;
        accettazioneOutput.job = accettazione.index;
        accettazioneOutput.lambda = (accettazione.index / (accettazione.lastArrival - accettazione.init_time));
        accettazioneOutput.Ts = accettazione.area / accettazione.index;
        accettazioneOutput.Tq = accettazione.queue / accettazione.index;
        accettazioneOutput.Es =  accettazione.service/ accettazione.index;
        accettazioneOutput.N = (accettazione.area) / (accettazione.current - accettazione.init_time);
        accettazioneOutput.Nq = accettazione.queue / (accettazione.current - accettazione.init_time);
        accettazioneOutput.ro = accettazione.service /(accettazione.current - accettazione.init_time);

        dati anamnesiOutput;
        anamnesiOutput.name = anamnesi.name;
        anamnesiOutput.job = anamnesi.index;
        anamnesiOutput.lambda = anamnesi.index / (anamnesi.lastArrival - anamnesi.init_time);
        anamnesiOutput.Ts = (anamnesi.area) / anamnesi.index;
        anamnesiOutput.Tq = anamnesi.queue / anamnesi.index;
        anamnesiOutput.Es = (anamnesi.service/anamnesi.servers) / anamnesi.index;
        anamnesiOutput.N = anamnesi.area / (anamnesi.current - anamnesi.init_time);
        anamnesiOutput.Nq = anamnesi.queue / (anamnesi.current - anamnesi.init_time);
        anamnesiOutput.ro = anamnesi.service / (anamnesi.servers*(anamnesi.current - anamnesi.init_time));

        dati inoculazioneOutput;
        inoculazioneOutput.name = inoculazione.name;
        inoculazioneOutput.job = inoculazione.index;
        inoculazioneOutput.lambda = inoculazione.index / (inoculazione.lastArrival - inoculazione.init_time);
        inoculazioneOutput.Ts = inoculazione.area / inoculazione.index;
        inoculazioneOutput.Tq = inoculazione.queue / inoculazione.index;
        inoculazioneOutput.Es = (inoculazione.service/inoculazione.servers) / inoculazione.index;
        inoculazioneOutput.N = (inoculazione.area) / (inoculazione.current - inoculazione.init_time);
        inoculazioneOutput.Nq = inoculazione.queue / (inoculazione.current- inoculazione.init_time);
        inoculazioneOutput.ro = inoculazione.service / (inoculazione.servers*(inoculazione.current - inoculazione.init_time));
        
        dati inoculazione_attenuatoOutput;
        inoculazione_attenuatoOutput.name = inoculazione_attenuato.name;
        inoculazione_attenuatoOutput.job = inoculazione_attenuato.index;
        inoculazione_attenuatoOutput.lambda = inoculazione_attenuato.index / (inoculazione_attenuato.lastArrival - inoculazione_attenuato.init_time);
        inoculazione_attenuatoOutput.Ts = inoculazione_attenuato.area / inoculazione_attenuato.index;
        inoculazione_attenuatoOutput.Tq = inoculazione_attenuato.queue / inoculazione_attenuato.index;
        inoculazione_attenuatoOutput.Es = (inoculazione_attenuato.service/inoculazione_attenuato.servers) / inoculazione_attenuato.index;
        inoculazione_attenuatoOutput.N = (inoculazione_attenuato.area) / (inoculazione_attenuato.current - inoculazione_attenuato.init_time);
        inoculazione_attenuatoOutput.Nq = inoculazione_attenuato.queue / (inoculazione_attenuato.current- inoculazione_attenuato.init_time);
        inoculazione_attenuatoOutput.ro = inoculazione_attenuato.service / (inoculazione_attenuato.servers*(inoculazione_attenuato.current - inoculazione_attenuato.init_time));
        datiFinite[0] = accettazioneOutput;
        datiFinite[1] = anamnesiOutput;
        datiFinite[2] = inoculazioneOutput;
        datiFinite[3] = inoculazione_attenuatoOutput;
        
    
    
    }
    
    return 0.0;
        
    }



double autocorr(double data[], size_t size) {
    
    double mean = 0;
    double variance = 0;
    double autocorr = 0;
    
    int i;

    
    for (i = 0; i < size; i++) {
        
        mean += data[i];
    }
    mean /= size;

    
    for (i = 0; i < size; i++) {
        variance += (data[i] - mean) * (data[i] - mean);
    }
    variance /= size;

 //for(int lag = 0; lag < size; lag ++){  
    
    for (i = 0; i < size - LAG; i++) {
        autocorr += (data[i] - mean) * (data[i + LAG] - mean);
    }
    autocorr/= (size - LAG) * variance;//}
    
    

    printf("L'autocorrelazione è : %.1f\n", autocorr);

    return 0;
}




double estimate(double statistics[], int size){

  long   n    = 0;                   
  double sum  = 0.0;
  double mean = 0.0;
  double data;
  double stdev;
  double u, t, w;
  double diff;
  

  for(int i=0; i<size; i++) {                 /* use Welford's one-pass method */
    data = statistics[i];                            /* to calculate the sample mean  */
    n++;                                        /* and standard deviation        */
    diff  = data - mean;
    sum  += diff * diff * (n - 1.0) / n ;
    mean += diff / n;
  }
  stdev  = sqrt(sum / n);

  if (n > 1) {
    u = 1.0 - (1.0 - LOC) / 2.0;              /* interval parameter  */
    t = idfStudent(n - 1, u);                 /* critical value of t */
    w = t * stdev / sqrt(n-1);              /* interval half width */
    printf("%.6f +/- %.6f\n", mean, w);
    
    }
    else printf("ERROR - insufficient data\n");
  return 0;
}
                              





int main(void){

    
    int digits;
    dati matrix[K][4]; 
                         
  
//accettazione 
    double lambda = LAMBDA;
    double Es = (1/MU_ACCETTAZIONE);
    double ro = lambda * Es;
    double Etq = (ro * Es) / (1-ro);
    double Ets = Etq + Es;


//anamnesi 
    double lambda_a = LAMBDA;  
    double Es_a = (1/(MU_ANAMNESI*SERVERS_ANAMNESI));
    double ro_a = lambda_a * Es_a;
    //double p_zero_accettazione_a = 1 / (  ((pow((SERVERS_ANAMNESI*ro_a),SERVERS_ANAMNESI)) / ((1-ro_a)*(fattoriale(SERVERS_ANAMNESI))) ) + 1 + (ro_a*SERVERS_ANAMNESI) + ((pow((SERVERS_ANAMNESI*ro_a),2)) / ((fattoriale(2))) + (ro_a*SERVERS_ANAMNESI) + ((pow((SERVERS_ANAMNESI*ro_a),3)) / ((fattoriale(3))) )));

    double p_zero_accettazione_a = 1 / (  ((pow((SERVERS_ANAMNESI*ro_a),SERVERS_ANAMNESI)) / ((1-ro_a)*(fattoriale(SERVERS_ANAMNESI))) ) + 1 + (ro_a*SERVERS_ANAMNESI) + ((pow((SERVERS_ANAMNESI*ro_a),2)) / ((fattoriale(2))) + (ro_a*SERVERS_ANAMNESI) + ((pow((SERVERS_ANAMNESI*ro_a),3)) / ((fattoriale(3))) )));
    //double p_zero_accettazione_a = 1 / (  ((pow((SERVERS_ANAMNESI*ro_a),SERVERS_ANAMNESI)) / ((1-ro_a)*(fattoriale(SERVERS_ANAMNESI))) ) + 1 + (ro_a*SERVERS_ANAMNESI) );
    double pq_accettazione_a = p_zero_accettazione_a * ((pow((SERVERS_ANAMNESI*ro_a),SERVERS_ANAMNESI)) / ((fattoriale(SERVERS_ANAMNESI))*(1-ro_a)));
    double Etq_a = ((pq_accettazione_a * (Es_a)) / (1 - ro_a));
    double Ets_a = Etq_a + (Es_a*SERVERS_ANAMNESI);
    


//inoculazione 
    double lambda_i = LAMBDA*(prob_no_vaccino)*(1-prob_vaccino_attenuato);
    double Es_i = (1/(MU_INOCULAZIONE*SERVERS_INOCULAZIONE));
    double ro_i = lambda_i * Es_i;
    //double p_zero_inoculazione = 1 / (  ((pow((SERVERS_INOCULAZIONE*ro_i),SERVERS_INOCULAZIONE)) / ((1-ro_i)*(fattoriale(SERVERS_INOCULAZIONE))) ) + 1 + (ro_i*SERVERS_INOCULAZIONE) + ((pow((SERVERS_INOCULAZIONE*ro_i),2)) / ((fattoriale(2))) ));
    double p_zero_inoculazione = 1 / (  ((pow((SERVERS_INOCULAZIONE*ro_i),SERVERS_INOCULAZIONE)) / ((1-ro_i)*(fattoriale(SERVERS_INOCULAZIONE))) ) + 1 + (ro_i*SERVERS_INOCULAZIONE) );
    double pq_inoculazione = p_zero_inoculazione * ((pow((SERVERS_INOCULAZIONE*ro_i),SERVERS_INOCULAZIONE)) / (fattoriale(SERVERS_INOCULAZIONE)*(1-ro_i)) );
    double Etq_i = ((pq_inoculazione * (Es_i)) / (1 - ro_i));
    double Ets_i = (Etq_i) + (Es_i*SERVERS_INOCULAZIONE);

 //inoculazione attenuato 
    double lambda_att = (LAMBDA*(prob_no_vaccino))*prob_vaccino_attenuato;
    double Es_att = (1/(MU_INOCULAZIONE_ATTENUATO*SERVERS_INOCULAZIONE_ATTENUATO));
    double ro_att = lambda_att * Es_att;
    double p_zero_inoculazione_att = 1 / (  ((pow((SERVERS_INOCULAZIONE_ATTENUATO*ro_att),SERVERS_INOCULAZIONE_ATTENUATO)) / ((1-ro_att)*(fattoriale(SERVERS_INOCULAZIONE_ATTENUATO))) ) + 1 + (ro_att*SERVERS_INOCULAZIONE_ATTENUATO));
    double pq_inoculazione_att = p_zero_inoculazione_att * ((pow((SERVERS_INOCULAZIONE_ATTENUATO*ro_att),SERVERS_INOCULAZIONE_ATTENUATO)) / (fattoriale(SERVERS_INOCULAZIONE_ATTENUATO)*(1-ro_att)) );
    double Etq_att = ((pq_inoculazione_att * (Es_att)) / (1 - ro_att));
    double Ets_att = (Etq_att) + (Es_att*SERVERS_INOCULAZIONE_ATTENUATO);
  
//


    printf("\n*******************************\n");
    printf("*                             *\n");
    printf("*       PROGETTO PMCSN        *\n");
    printf("*      CENTRO VACCINALE       *\n");
    printf("*                             *\n");
    printf("*******************************\n");

    PlantSeeds(SEED);
    
    
    printf("\nSelezionare(1/2):\n");
    printf("\n\033[33m1)SIMULAZIONE\033[0m\n");
    printf("\n\033[33m2)ANALISI ({1,2,2,2})\033[0m\n");

    if(SERVERS_ACCETTAZIONE > 1){
        
            printf("\033[33mIl centro accettazione è un centro SINGLE, diminuire il numero di server\033[0m\n");
            system( "read -n 1 -s -p \"Premere un tasto per continuare...\"" );
 
            return 0;
            
            
        }

    if(SERVERS_ANAMNESI == 1){
        
            printf("\033[33mIl centro anamnesi è un centro MULTI, aumentare il numero di server\033[0m\n");
            system( "read -n 1 -s -p \"Premere un tasto per continuare...\"" );
 
            return 0;
        }
    
    if(SERVERS_INOCULAZIONE == 1){
        
            printf("\033[33mIl centro inoculazione è un centro MULTI, aumentare il numero di server\033[0m\n");
            system( "read -n 1 -s -p \"Premere un tasto per continuare...\"" );
 
            return 0;
        }
    
    if(SERVERS_INOCULAZIONE_ATTENUATO == 1){
        
            printf("\033[33mIl centro inoculazione_attenuato è un centro MULTI, aumentare il numero di server\033[0m\n");
            system( "read -n 1 -s -p \"Premere un tasto per continuare...\"" );
 
            return 0;
        }
    
    scanf("%d",&digits);
    if(digits == 1){
        
        printf("\n\033[33m1)ORIZZONTE FINITO\033[0m\n");
        printf("\n\033[33m2)ORIZZONTE INFINITO\033[0m\n");
        scanf("%d",&digits);
        if(digits == 1){
        
        printf("\033[33m\nSIMULAZIONE AD ORIZZONTE FINITO\n\033[0m\n");

           
        printf("-----------------------------SIMULAZIONE----------------------------- \n\n");
 

            for(int i=0; i < ITERATIONS; i++){
            printf(" ");
            sim("finite", matrix[i], NULL);
            
        }
//accettazione
        long n[7];
        double mean[7];
        double sum[7];
        double data[7];
        double stdev[7];
        double t[7];
        double w[7];
        double u[7];
        double diff[7];
//attenuato
        double diffatt[7];
        double meanatt[7];
        double sumatt[7];
        double stdevatt[7];
        double natt[7];
        double tatt[7];
        double watt[7];
//anamnesi
        double diffan[7];
        double meanan[7];
        double suman[7];
        double stdevan[7];
        double nana[7];
        double tan[7];
        double wan[7];
//inoculazione
        double diffin[7];
        double meanin[7];
        double sumin[7];
        double stdevin[7];
        double nin[7];
        double tin[7];
        double win[7];

        printf("\n\n\n\033[33mCENTRO DI ACCETTAZIONE \n\n\033[0m " );

                
                for(int i=0; i < ITERATIONS - 1; i++) { 
                    
                                        
                    for(int z = 0; z<7; z++){   
                        n[z]++;     }                            
                        diff[0]  = matrix[i][0].lambda - mean[0];
                        diff[1]  = matrix[i][0].Ts - mean[1];
                        diff[2] = matrix[i][0].Tq - mean[2];
                        diff[3] = matrix[i][0].Es - mean[3];
                        diff[4] = matrix[i][0].N - mean[4];
                        diff[5] = matrix[i][0].Nq - mean[5];
                        diff[6] = matrix[i][0].ro - mean[6];
                        

                    for(int i=0; i<7; i++){
                        sum[i]  += diff[i] * diff[i] * (n[i] - 1.0) / n[i];
                        mean[i] += diff[i] / n[i];
                        stdev[i]  = sqrt(sum[i] / n[i]);
                                                    
                        t[i] = idfStudent(n[i] - 1, 1.0 - (1.0 - LOC) / 2.0);                      
                        w[i] = t[i] * stdev[i] / sqrt(n[i] - 1);

                    }
                } 

                printf("************************************************************\n");
                        printf("*  λ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[0], w[0]);
                        printf("*  E(Ts)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[1], w[1]);
                        printf("*  E(Tq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[2], w[2]);
                        printf("*  E(S)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[3], w[3]);
                        printf("*  E(N)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[4], w[4]);
                        printf("*  E(Nq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[5], w[5]);
                        printf("*  ρ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[6], w[6]);
                        printf("************************************************************\n");
                        printf("\nVerifica di consistenza (E(Ts) = E(Tq) + E(s)): %f\n", (mean[2] + (mean[3]*SERVERS_ACCETTAZIONE) - mean[1]));

                printf("\n\n\n\033[33mCENTRO DI ANAMNESI \n\n\033[0m " );
                for(int i=0; i < ITERATIONS - 1; i++) { 
                    
                                        
                    for(int z = 0; z<7; z++){   
                        nana[z]++;     }                            
                        diffan[0]  = matrix[i][1].lambda - meanan[0];
                        diffan[1]  = matrix[i][1].Ts - meanan[1];
                        diffan[2] = matrix[i][1].Tq - meanan[2];
                        diffan[3] = matrix[i][1].Es - meanan[3];
                        diffan[4] = matrix[i][1].N - meanan[4];
                        diffan[5] = matrix[i][1].Nq - meanan[5];
                        diffan[6] = matrix[i][1].ro - meanan[6];
                        

                    for(int i=0; i<7; i++){
                        suman[i]  += diffan[i] * diffan[i] * (nana[i] - 1.0) / nana[i];
                        meanan[i] += diffan[i] / nana[i];
                        stdevan[i]  = sqrt(suman[i] / nana[i]);
                                                    
                        tan[i] = idfStudent(nana[i] - 1, 1.0 - (1.0 - LOC) / 2.0);                      
                        wan[i] = tan[i] * stdevan[i] / sqrt(nana[i] - 1);

                    }
                } 

                printf("************************************************************\n");
                        printf("*  λ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[0], wan[0]);
                        printf("*  E(Ts)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[1], wan[1]);
                        printf("*  E(Tq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[2], wan[2]);
                        printf("*  E(S)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[3], wan[3]);
                        printf("*  E(N)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[4], wan[4]);
                        printf("*  E(Nq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[5], wan[5]);
                        printf("*  ρ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[6], wan[6]);
                        printf("************************************************************\n");
                        printf("\nVerifica di consistenza (E(Ts) = E(Tq) + E(s)): %f\n", (meanan[2] + (meanan[3]*SERVERS_ANAMNESI) - meanan[1]));

                printf("\n\n\n\033[33mCENTRO DI INOCULAZIONE\n\n\033[0m " );
                for(int i=0; i < ITERATIONS - 1; i++) { 
                    
                                        
                    for(int z = 0; z<7; z++){   
                        nin[z]++;     }                            
                        diffin[0]  = matrix[i][2].lambda - meanin[0];
                        diffin[1]  = matrix[i][2].Ts - meanin[1];
                        diffin[2] = matrix[i][2].Tq - meanin[2];
                        diffin[3] = matrix[i][2].Es - meanin[3];
                        diffin[4] = matrix[i][2].N - meanin[4];
                        diffin[5] = matrix[i][2].Nq - meanin[5];
                        diffin[6] = matrix[i][2].ro - meanin[6];
                        

                    for(int i=0; i<7; i++){
                        sumin[i]  += diffin[i] * diffin[i] * (nin[i] - 1.0) / nin[i];
                        meanin[i] += diffin[i] / nin[i];
                        stdevin[i]  = sqrt(sumin[i] / nin[i]);
                                                    
                        tin[i] = idfStudent(nin[i] - 1, 1.0 - (1.0 - LOC) / 2.0);                      
                        win[i] = tin[i] * stdevin[i] / sqrt(nin[i] - 1);
 
                    }
                } 

                printf("************************************************************\n");
                        printf("*  λ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[0], win[0]);
                        printf("*  E(Ts)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[1], win[1]);
                        printf("*  E(Tq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[2], win[2]);
                        printf("*  E(S)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[3], win[3]);
                        printf("*  E(N)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[4], win[4]);
                        printf("*  E(Nq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[5], win[5]);
                        printf("*  ρ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[6], win[6]);
                        printf("************************************************************\n");
                        printf("\nVerifica di consistenza (E(Ts) = E(Tq) + E(s)): %f\n", (meanin[2] + (meanin[3]*SERVERS_INOCULAZIONE) - meanin[1]));
                
                printf("\n\n\n\033[33mCENTRO DI INOCULAZIONE ATTENUATO \n\n\033[0m " );
                for(int i=0; i < ITERATIONS - 1; i++) { 
                

                                        
                    for(int z = 0; z<7; z++){   
                        natt[z]++;     }                            
                        diffatt[0]  = matrix[i][3].lambda - meanatt[0];
                        diffatt[1]  = matrix[i][3].Ts - meanatt[1];
                        diffatt[2] = matrix[i][3].Tq - meanatt[2];
                        diffatt[3] = matrix[i][3].Es - meanatt[3];
                        diffatt[4] = matrix[i][3].N - meanatt[4];
                        diffatt[5] = matrix[i][3].Nq - meanatt[5];
                        diffatt[6] = matrix[i][3].ro - meanatt[6];
                        

                    for(int i=0; i<7; i++){
                        sumatt[i]  += diffatt[i] * diffatt[i] * (natt[i] - 1.0) / natt[i];
                        meanatt[i] += diffatt[i] / natt[i];
                        stdevatt[i]  = sqrt(sumatt[i] / natt[i]);
                                                    
                        tatt[i] = idfStudent(natt[i] - 1, 1.0 - (1.0 - LOC) / 2.0);                      
                        watt[i] = tatt[i] * stdevatt[i] / sqrt(natt[i] - 1);
                        
 
                    }
                } 

                printf("************************************************************\n");
                        printf("*  λ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[0], watt[0]);
                        printf("*  E(Ts)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[1], watt[1]);
                        printf("*  E(Tq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[2], watt[2]);
                        printf("*  E(S)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[3], watt[3]);
                        printf("*  E(N)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[4], watt[4]);
                        printf("*  E(Nq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[5], watt[5]);
                        printf("*  ρ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[6], watt[6]);
                        printf("************************************************************\n");
                        printf("\nVerifica di consistenza (E(Ts) = E(Tq) + E(s)): %f\n", (meanatt[2] + (meanatt[3]*SERVERS_INOCULAZIONE_ATTENUATO) - meanatt[1]));

            
        
        double Ts_finite[ITERATIONS];
        double Tq_finite[ITERATIONS];
        for(int i=0; i<ITERATIONS; i++){

            Ts_finite[i] = (matrix[i][0].Ts)+
                matrix[i][1].Ts + (prob_no_vaccino)*(matrix[i][2].Ts)*(1-prob_vaccino_attenuato) + (prob_no_vaccino)*(matrix[i][3].Ts)*(prob_vaccino_attenuato);
            Tq_finite[i] = (matrix[i][0].Tq)+
                matrix[i][1].Tq + (prob_no_vaccino)*(matrix[i][2].Tq)*(1-prob_vaccino_attenuato) + (prob_no_vaccino)*(matrix[i][3].Tq)*(prob_vaccino_attenuato);                   
        }
        
        printf("\n\n");
        printf("\n\n\033[33mProbabilità di utenti che effettuano il vaccino:  \n\n\033[0m %.6f %% ", prob_no_vaccino);
        printf("\n\n\033[33mProbabilità di pazienti che fanno vaccino vivo:  \n\n\033[0m %.6f %%", prob_vaccino_attenuato);
        printf("\n\n\033[33mTempo di attesa del centro vaccinale: \n\033[0m\n");
        estimate(Tq_finite, ITERATIONS);
        printf("\n\n\033[33mTempo di risposta del centro vaccinale: \n\033[0m\n");
        estimate(Ts_finite, ITERATIONS);
        printf("\n");

        }
    else if(digits == 2){
        
        
        printf("\033[33m\nSIMULAZIONE AD ORIZZONTE INFINITO\n\033[0m\n"); 
        printf("-----------------------------SIMULAZIONE----------------------------- \n\n");
        sim("infinite", NULL, matrix);
// accettazione
        long n[7];  
        double mean[7];
        double sum[7];
        double data[7];
        double stdev[7];
        double t[7];
        double w[7];
        double u[7];
        double diff[7];
//inoculazione attenuato
        double diffatt[7];
        double meanatt[7];
        double sumatt[7];
        double stdevatt[7];
        double natt[7];
        double tatt[7];
        double watt[7];
//anamnesi
        double diffan[7];
        double meanan[7];
        double suman[7];
        double stdevan[7];
        double nana[7];
        double tan[7];
        double wan[7];
//inoculazione
        double diffin[7];
        double meanin[7];
        double sumin[7];
        double stdevin[7];
        double nin[7];
        double tin[7];
        double win[7];
        
        printf("\n\n\n\033[33mCENTRO DI ACCETTAZIONE \n\n\033[0m " );
        
        
                
            for(int i=0; i < K; i++) { 
                     
                    diff[0]  = matrix[i][0].lambda - mean[0];
                    diff[1]  = matrix[i][0].Ts - mean[1];
                    diff[2] = matrix[i][0].Tq - mean[2];
                    diff[3] = matrix[i][0].Es - mean[3];
                    diff[4] = matrix[i][0].N - mean[4];
                    diff[5] = matrix[i][0].Nq - mean[5];
                    diff[6] = matrix[i][0].ro - mean[6];
                        

                    for(int i=0; i<7; i++){
                        n[i]++;
                        sum[i]  += diff[i] * diff[i] * (n[i] - 1.0) / n[i];
                        mean[i] += diff[i] / n[i];
                        stdev[i]  = sqrt(sum[i] / n[i]);
                        if(n[i]>1){
                                                    
                        t[i] = idfStudent(n[i] - 1, 1.0 - (1.0 - LOC) / 2.0);                      
                        w[i] = t[i] * stdev[i] / sqrt(n[i] - 1);}
                       
                        

                    }
                } 

                
 printf(" based upon %ld data points", n[0]);
                printf("************************************************************\n");
                        printf("*  λ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[0], w[0]);
                        printf("*  E(Ts)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[1], w[1]);
                        printf("*  E(Tq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[2], w[2]);
                        printf("*  E(S)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[3], w[3]);
                        printf("*  E(N)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[4], w[4]);
                        printf("*  E(Nq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[5], w[5]);
                        printf("*  ρ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", mean[6], w[6]);
                        printf("************************************************************\n");
                        printf("\nVerifica di consistenza (E(Ts) = E(Tq) + E(s)): %f\n", (mean[2] + (mean[3]*SERVERS_ACCETTAZIONE) - mean[1]));

                printf("\n\n\n\033[33mCENTRO DI ANAMNESI \n\n\033[0m " );
                for(int i=0; i < K; i++) { 
                    
                                        
                    for(int z = 0; z<7; z++){   
                        nana[z]++;     }                            
                        diffan[0]  = matrix[i][1].lambda - meanan[0];
                        diffan[1]  = matrix[i][1].Ts - meanan[1];
                        diffan[2] = matrix[i][1].Tq - meanan[2];
                        diffan[3] = matrix[i][1].Es - meanan[3];
                        diffan[4] = matrix[i][1].N - meanan[4];
                        diffan[5] = matrix[i][1].Nq - meanan[5];
                        diffan[6] = matrix[i][1].ro - meanan[6];
                        

                    for(int i=0; i<7; i++){
                        suman[i]  += diffan[i] * diffan[i] * (nana[i] - 1.0) / nana[i];
                        meanan[i] += diffan[i] / nana[i];
                        stdevan[i]  = sqrt(suman[i] / nana[i]);
                                                    
                        tan[i] = idfStudent(nana[i] - 1, 1.0 - (1.0 - LOC) / 2.0);                      
                        wan[i] = tan[i] * stdevan[i] / sqrt(nana[i] - 1);
                        

                    }
                } 
                

                printf("************************************************************\n");
                        printf("*  λ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[0], wan[0]);
                        printf("*  E(Ts)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[1], wan[1]);
                        printf("*  E(Tq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[2], wan[2]);
                        printf("*  E(S)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[3], wan[3]);
                        printf("*  E(N)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[4], wan[4]);
                        printf("*  E(Nq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[5], wan[5]);
                        printf("*  ρ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanan[6], wan[6]);
                        printf("************************************************************\n");
                        printf("\nVerifica di consistenza (E(Ts) = E(Tq) + E(s)): %f\n", (meanan[2] + (meanan[3]*SERVERS_ANAMNESI) - meanan[1]));

                printf("\n\n\n\033[33mCENTRO DI INOCULAZIONE\n\n\033[0m " );
                for(int i=0; i < 2; i++) { 
                    
                                        
                    for(int z = 0; z<7; z++){   
                        nin[z]++;     }                            
                        diffin[0]  = matrix[i][2].lambda - meanin[0];
                        diffin[1]  = matrix[i][2].Ts - meanin[1];
                        diffin[2] = matrix[i][2].Tq - meanin[2];
                        diffin[3] = matrix[i][2].Es - meanin[3];
                        diffin[4] = matrix[i][2].N - meanin[4];
                        diffin[5] = matrix[i][2].Nq - meanin[5];
                        diffin[6] = matrix[i][2].ro - meanin[6];
                        

                    for(int i=0; i<K; i++){
                        sumin[i]  += diffin[i] * diffin[i] * (nin[i] - 1.0) / nin[i];
                        meanin[i] += diffin[i] / nin[i];
                        stdevin[i]  = sqrt(sumin[i] / nin[i]);
                                                    
                        tin[i] = idfStudent(nin[i] - 1, 1.0 - (1.0 - LOC) / 2.0);                      
                        win[i] = tin[i] * stdevin[i] / sqrt(nin[i] - 1);

                    }
                } 

                printf("************************************************************\n");
                        printf("*  λ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[0], win[0]);
                        printf("*  E(Ts)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[1], win[1]);
                        printf("*  E(Tq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[2], win[2]);
                        printf("*  E(S)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[3], win[3]);
                        printf("*  E(N)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[4], win[4]);
                        printf("*  E(Nq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[5], win[5]);
                        printf("*  ρ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanin[6], win[6]);
                        printf("************************************************************\n");
                        printf("\nVerifica di consistenza (E(Ts) = E(Tq) + E(s)): %f\n", (meanin[2] + (meanin[3]*SERVERS_INOCULAZIONE) - meanin[1]));
                
                printf("\n\n\n\033[33mCENTRO DI INOCULAZIONE ATTENUATO \n\n\033[0m " );
                for(int i=0; i < K; i++) { 
                

                                        
                    for(int z = 0; z<7; z++){   
                        natt[z]++;     }                            
                        diffatt[0]  = matrix[i][3].lambda - meanatt[0];
                        diffatt[1]  = matrix[i][3].Ts - meanatt[1];
                        diffatt[2] = matrix[i][3].Tq - meanatt[2];
                        diffatt[3] = matrix[i][3].Es - meanatt[3];
                        diffatt[4] = matrix[i][3].N - meanatt[4];
                        diffatt[5] = matrix[i][3].Nq - meanatt[5];
                        diffatt[6] = matrix[i][3].ro - meanatt[6];
                        

                    for(int i=0; i<7; i++){
                        sumatt[i]  += diffatt[i] * diffatt[i] * (natt[i] - 1.0) / natt[i];
                        meanatt[i] += diffatt[i] / natt[i];
                        stdevatt[i]  = sqrt(sumatt[i] / natt[i]);
                                                    
                        tatt[i] = idfStudent(natt[i] - 1, 1.0 - (1.0 - LOC) / 2.0);                      
                        watt[i] = tatt[i] * stdevatt[i] / sqrt(natt[i] - 1);

                    }
                } 

                printf("************************************************************\n");
                        printf("*  λ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[0], watt[0]);
                        printf("*  E(Ts)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[1], watt[1]);
                        printf("*  E(Tq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[2], watt[2]);
                        printf("*  E(S)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[3], watt[3]);
                        printf("*  E(N)                      =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[4], watt[4]);
                        printf("*  E(Nq)                     =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[5], watt[5]);
                        printf("*  ρ                         =   ");
                        printf("%.6f +/- %.6f         \n\n", meanatt[6], watt[6]);
                        printf("************************************************************\n");
                        printf("\nVerifica di consistenza (E(Ts) = E(Tq) + E(s)): %f\n", (meanatt[2] + (meanatt[3]*SERVERS_INOCULAZIONE_ATTENUATO) - meanatt[1]));

           
        
        
        double Ts_infinite[K];
        double Tq_infinite[K];
        for(int i=0; i<K; i++){
            
                        Ts_infinite[i] = (matrix[i][0].Ts)+
                            matrix[i][1].Ts + (prob_no_vaccino)*(matrix[i][2].Ts)*(1-prob_vaccino_attenuato) + (prob_no_vaccino)*(matrix[i][3].Ts)*(prob_vaccino_attenuato);
                        Tq_infinite[i] = (matrix[i][0].Tq)+
                            matrix[i][1].Tq + (prob_no_vaccino)*(matrix[i][2].Tq)*(1-prob_vaccino_attenuato) + (prob_no_vaccino)*(matrix[i][3].Tq)*(prob_vaccino_attenuato);
            
            
        }
        printf("\n\n\n\033[33mAutocorrelazione dei dati nel campione utilizzando K = %d e B = %d: \n\n\033[0m ", K,B );
        autocorr(Ts_infinite,K);
        printf("\n\n\033[33mProbabilità di utenti che effettuano il vaccino:  \n\n\033[0m %.6f %%", prob_no_vaccino);
        printf("\n\n\033[33mProbabilità di pazienti che fanno vaccino vivo:  \n\n\033[0m %.6f %%", prob_vaccino_attenuato);
        printf("\n\n\n\033[33mTempo di attesa del centro vaccinale: \n\n\033[0m " );
        estimate(Tq_infinite, K);
        printf("\n\n\n\033[33mTempo di risposta del centro vaccinale: \n\n\033[0m " );
        estimate(Ts_infinite, K);
        printf("\n\n");
    }
    }else if(digits == 2){
        printf("-------------------------------ANALISI------------------------------- \n\n");
            

            printf("\nACCETTAZIONE\n\n");
            printf("**************************************************\n");
            printf("*  λ                         = %6.6f req/min  *\n", LAMBDA);            
            printf("*  E(S)                      = %6.6f min      *\n", (1/(MU_ACCETTAZIONE)));
            printf("*  ρ                         = %6.6f min      *\n", (lambda * Es));
            printf("*  E(Tq)                     = %6.6f min      *\n", Etq);
            printf("*  E(Ts)                     = %6.6f min      *\n", Etq + (Es));           
            printf("**************************************************\n");

           

            printf("\nANAMNESI\n\n");
            printf("**************************************************\n");
            printf("*  λ                         = %6.6f req/min  *\n", lambda_a);            
            printf("*  E(S)                      = %6.6f min      *\n", Es_a);
            printf("*  ρ                         = %6.6f min      *\n", (lambda_a * Es_a));
            printf("*  E(Tq)                     = %6.6f min      *\n", Etq_a);
            printf("*  E(Ts)                     = %6.6f min      *\n", Ets_a);           
            printf("**************************************************\n");
            

            printf("\nINOCULAZIONE\n\n");
            printf("**************************************************\n");
            printf("*  λ                         = %6.6f req/min  *\n", lambda_i);             
            printf("*  E(S)                      = %6.6f min      *\n", Es_i);
            printf("*  ρ                         = %6.6f min      *\n", lambda_i * Es_i);
            printf("*  E(Tq)                     = %6.6f min      *\n", Etq_i);
            printf("*  E(Ts)                     = %6.6f min      *\n", Ets_i);           
            printf("**************************************************\n"); 

            printf("\nINOCULAZIONE ATTENUATO\n\n");
            printf("**************************************************\n");
            printf("*  λ                         = %6.6f req/min  *\n", lambda_att);             
            printf("*  E(S)                      = %6.6f min      *\n", Es_att);
            printf("*  ρ                         = %6.6f min      *\n", ro_att);
            printf("*  E(Tq)                     = %6.6f min      *\n", Etq_att);
            printf("*  E(Ts)                     = %6.6f min      *\n", Ets_att);           
            printf("**************************************************\n");
            printf("\n\n");
            printf("\033[33mTempo di attesa del centro vaccinale: %6.6f\n\n\n\033[0m", (Etq) + Etq_a + (prob_no_vaccino)*Etq_i*(1-prob_vaccino_attenuato) +(prob_no_vaccino)*Etq_att*(prob_vaccino_attenuato) );

            printf("\n\n\033[33mTempo di risposta del centro vaccinale: %6.6f\n\n\n\033[0m", (Ets) + Ets_a + (prob_no_vaccino)*Ets_i*(1-prob_vaccino_attenuato) +(prob_no_vaccino)*Ets_att*(prob_vaccino_attenuato) );
    }

    return 0;
}

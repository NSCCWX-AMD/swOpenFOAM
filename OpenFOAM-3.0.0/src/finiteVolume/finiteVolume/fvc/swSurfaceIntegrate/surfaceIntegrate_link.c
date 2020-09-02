//link host and slave
    #include <athread.h>
    #include"surfaceIntegrate_struct.h"
    #include<stdio.h>
     extern void SLAVE_FUN(surfaceIntegrate_slave)(struct surfaceIntegrate_para*);

//static if_owner_init=1;
void SWlink_slave(struct surfaceIntegrate_para *para){
   /*  if(if_owner_init){
     athread_init();
	 athread_enter64();
     if_owner_init=0;
     }*/
     //athread_spawn64(slave_surfaceIntegrate_slave,para);
     //athread_join64();
     athread_spawn(surfaceIntegrate_slave,para);
     athread_join();

}















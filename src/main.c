# include "tsunami.h"

int main(void)
{   
        char *meshName = "../data/PacificFine.txt";
        char *resultBaseName = "output/tsunamiFine";
        //tsunamiCompute(2,30000,50,meshName,resultBaseName); //T'augmente le pas de temps si tu prends un maillage plus grossier  2 pour fine 20 pour Tiny
        tsunamiAnimate(2,30000,500,meshName,resultBaseName);
        exit(0);     
}

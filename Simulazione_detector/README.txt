File istruzioni per simulazione montecarlo Zambon-Meinardi

.L Compile.C, dovrebbe compilare automaticamente tutte e 5 le classi nell'ordine corretto. (TracciaMC.cxx, Vertex.cxx, Multis.cxx, Smearing.cxx, Punto2.cxx)
per runnare la simulazione bisogna prima loadare il file con .L Simulazione2.C 
e poi chiamare la funzione Simulazione2(b1,b2,b3), dove gli argomenti b1,b2,b3 sono booleani:
- b1 per il multiple scattering (0 = non attivo, 1 = attivo)
- b2 per la distribuzione di molteplicità (0 = uniforme, 1 = distribuzione data)
- b3 per la distribuzione di eta (0 = uniforme, 1 = distribuzione data)
al termine della simulazione, accendendo al TBrowser, si possono vedere le distribuzioni di molteplicità
e le coordinate di ciascun vertice simulato, salvati come foglie del TTree2 nel ramo "vertMult", e le coordinate delle hits sui vari layers,
salvate sugli altri rami.
Per eseguire la ricostruzione è sufficiente .x Ricostruzione2.C; al termine di essa sarà possibile visualizzare,
sempre sul TBrowser, l'istogramma dei residui salvato nel file "zrec-ztrue-tree.root".
Infine .x Eff.C esegue la macro relativa ad efficienza e risoluzione che utilizza i dati dei tree 3 e 4. Vengono eseguiti 
i fit gaussiani sugli istogrammi dei residui per un dato valore di molteplicità, da cui si estrae la sigma per la risoluzione.
I grafici di efficienza e risoluzione vengono salvati nel file "efficienza.root".
 







dubbi: quanti punti spuri? grafico efficienza e risol per entrambe le distrib di moltep?
// Your First C++ Program
using namespace std;
#include <iostream>
#include <vector>
#include <map>
#include <tuple>


double state_probability (int i, int j, int k);
double calcolo_lambda_j();

int i = 18;
int j = 20;
int k = 0;
const int Number_of_nodes = 150;
const double lambda = 2;
const int mu_j = 1;
const bool type_of_node = 0; //0 per benevolo, 1 per malevolo
const int max_allocab_resources = 3;
const int max_passo = j+7; //esiste una formuletta per calcolarlo

tuple <int, int, int> stato;
map<tuple<int, int, int>, double> mapOfTuple;

int main() {
   
  
    stato = make_tuple(i, j, k);
    
    double P = 1;
    mapOfTuple[stato] = P;
    cout << "Stato: (" << i << "," << j << "," << k << ") Il valore di Probabilita' di stato è: " << P << endl;
  //esempi
  //  stato = make_tuple(18, 21, 0);
   // mapOfTuple[stato] = 3;
 //    int estract = mapOfTuple[make_tuple(i,j,k)];


    //la prima prob è uno, quindi inizio a calcolare dalla seconda
    k++;
 //   stato = make_tuple(i, j, k);
 //   mapOfTuple[stato] = P;
 //   k++;

    for (j=20; j <= max_passo; j++) {
        for (i=18; i <= j-2; i++) {
            for (k=0; k <= max_allocab_resources; k++) {
                if (j == 20 && i == 18 && k == 0)
                    continue;
                cout  <<"Stato: (" << i << "," << j << "," << k << ") ";
                stato = make_tuple(i, j, k);
                P = state_probability(i, j, k);
                mapOfTuple[stato] = P;
                cout << "Il valore di Probabilita' di stato è: " << P << endl;
             
            }
        }
    }
  

    return 0;
}

double state_probability(int i, int j, int k) {
   
    double probabilita = 0;
    double probabilita_feedback_positivo = 0;
    double probabilita_feedback_negativo = 0;
    double lambda_j = calcolo_lambda_j();
   
    if (type_of_node == 0) {
        probabilita_feedback_positivo = 0.9;
    }
    else {
        probabilita_feedback_positivo = 0.5;
    }

    probabilita_feedback_negativo = (1 - probabilita_feedback_positivo);

    // al posto di 0.5 andrà rispettivamente il valore dal map P(kj,Tj,nj-1), P(kj,Tj-1,nj+1), P(kj-1,Tj-1,nj+1)
    //ricordasi di sostuire tutti e tre gli act_res[k] rispettivamente con k-1, e poi tutti k+1 quando si potrà fare

    //controllo sullo stato attraverso map

    double prob_prima_parte = 0;
    if (k - 1 >= 0){
        prob_prima_parte = mapOfTuple[make_tuple(i, j, k - 1)] * lambda_j / (lambda_j + (k - 1) * mu_j);
    }
    else{
        prob_prima_parte = 0;
    }
        
    double prob_seconda_parte = mapOfTuple[make_tuple(i, j-1, k+1)] * ((k + 1) * mu_j * probabilita_feedback_negativo) / (lambda_j + (k + 1) * mu_j);
    double prob_terza_parte = mapOfTuple[make_tuple(i-1, j-1, k+1)] * ((k + 1) * mu_j * probabilita_feedback_positivo) / (lambda_j + (k + 1) * mu_j);

    probabilita = prob_prima_parte + prob_seconda_parte + prob_terza_parte;
    return probabilita;
}
 

double calcolo_lambda_j() {
    double lambda_provider = 1;
    double prob_richiesta_assegnata_a_j = 1;
    double prob_soprasoglia = 1;
    double prob_j_piu_trusted = 1;
    double prob_j_non_piu_trusted_ma_risorse = 1;

   //prob_soprasoglia = ;
   // prob_j_piu_trusted = 1;
   // prob_j_non_piu_trusted_ma_risorse = 1;

    prob_richiesta_assegnata_a_j = prob_soprasoglia * prob_j_piu_trusted * prob_j_non_piu_trusted_ma_risorse ;
    lambda_provider = (lambda / Number_of_nodes)* prob_richiesta_assegnata_a_j;

    return lambda_provider;
}
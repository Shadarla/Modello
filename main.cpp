// Your First C++ Program
using namespace std;
#include <iostream>
#include <vector>
#include <map>
#include <tuple>


double state_probability (vector<int> pos_feed, vector<int> tot_feed, vector<int> act_res, int i, int j, int k);

const int Number_of_nodes = 150;
const int lambda = 2;
const int mu_j = 1;
const bool type_of_node = 0; //0 per benevolo, 1 per malevolo
const int max_allocab_resources = 3;
const int max_passo = 7; //esiste una formuletta per calcolarlo

tuple <int, int, int> stato;
map<tuple<int, int, int>, int> mapOfTuple;

int main() {
   
    stato = make_tuple(18, 20, 0);
    
    vector<int> positive_feedback, total_feedback, actual_resources;
    vector<double> Delta;
    double P = 0;
    mapOfTuple[stato] = P;

    int i = 20;
    int j = 20;
    int k = 2;
    stato = make_tuple(18, 21, 0);
    mapOfTuple[stato] = 3;
    stato = make_tuple(18, 20, 2);
    mapOfTuple[stato] = 4;
    stato = make_tuple(19, 22, 1);
    mapOfTuple[stato] = 1;
    int estract = mapOfTuple[make_tuple(i,j,k)];

    positive_feedback.push_back(18);
    total_feedback.push_back(20);
    actual_resources.push_back(0);

    //indici temporanei poi li sostituisco con i for

    for (j; j < max_passo; j++) {
        for (i; i <= j; i++) {
            for (k; k <= max_allocab_resources; k++) {
                P = state_probability(positive_feedback, total_feedback, actual_resources, i, j, k);
                cout << "Il valore di Probabilita' di stato al passo " << j << "e': " << P << endl;
                //actual_resources.push_back(k+1);

            }
            //positive_feedback.push_back(i);
        }
        //total_feedback.push_back(j);
    }
   
    //struct o map per lo stato da restituire alla funzione

    return 0;
}

double state_probability(vector<int> pos_feed, vector<int> tot_feed, vector<int> act_res, int i, int j, int k) {
   
    double probabilita = 0;
    double probabilita_feedback_positivo = 0;
    double probabilita_feedback_negativo = 0;
    int lambda_j = 2;

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

    probabilita = 0.5 * lambda_j / (lambda_j + (act_res[k-1] * mu_j)) + 0.5 * (act_res[k] * mu_j * probabilita_feedback_negativo) / (lambda_j + (act_res[k] * mu_j)) + 0.5 * (act_res[k] * mu_j * probabilita_feedback_positivo) / (lambda_j + (act_res[k] * mu_j));
    return probabilita;
}
// Modello
using namespace std;
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>


//C:\Users\gianc\Documents\GitHub\SSIoT\Sim - n_services_1 - n_devices_25 - n_master_1 - lambda_10.000000 - tot_sim_500 - seed_3 - resource_ctrl_1 - qoe_ctrl_1

double state_probability (int i, int j, int k);
double calcolo_lambda_j();
double trust_model(int id_amico_di_j);
void calcolo_denominatore_lambda_k();
double calcolo_lambda_k(int id_nodo_k);
double loss_probability();
double prob_blocco(int index);
double factorial(int n);
double prob_binomiale(vector<int> amici_di_i, int id_amico_j);
double prob_binomiale(vector<int> amici_di_i, int id_amico_j, int j);
int proporzionalita_tk(specifiche_nodo km);
double prob_di_blocco_generica(specifiche_nodo km, double lambda_k);

int kj = 18;
int Tj = 20;
int nj = 0;
const int Number_of_nodes = 25;
const double lambda = 20;
const int num_classi_di_servizio = 2;
const double mu_j = 1.4; //verifica (potrebbe anche non servire piu)
const bool type_of_node = 0; //0 per benevolo, 1 per malevolo (potrebbe anche non servire piu)
const int max_allocab_resources = 2;
const int max_passo = Tj+7; //esiste una formuletta per calcolarlo
const int num_amici = 13; // supposto uguale per tutti (to check per caso non omogeneo e potrebbe non servire piu)
double denominatore_lambda_k = 0; // non varia mai
const int id_nodo_j = 1; //id del nodo da valutare
const double lambda_ij = (lambda / Number_of_nodes); //essendo omogenea la distribuzione al momento è uguale per tutti
const double soglia = 0.03;


//vettori per prob blocco
vector<double> mu;
vector<double> num_amici_per_classe;
vector<double> array_prob_blocco_per_classe;
vector<double> Lambda;


class specifiche_nodo
{
private:
    int id_nodo;
    int num_amici;
    int classe; 
    double mu;
    bool type_of_node;
    double probabilita_feedback_positivo;
     
public:
    vector<double> S; //social factor
    specifiche_nodo();
   

    void set_id_nodo(int id) {
        id_nodo = id;
    }

    int get_id_nodo() {
        return this->id_nodo;
    }

    void set_num_amici(int amici_nodo) {
        num_amici = amici_nodo;
    }

    int get_num_amici() {
        return this->num_amici;
    }

    void set_type(bool tipo) {
        type_of_node = tipo;
    }

    int get_type() {
        return this->type_of_node;
    }

    void set_probabilità_feedback_positivo(bool type_of_node) {
        if (type_of_node == 0) {
            //nodo benevolo
            probabilita_feedback_positivo = 0.9;
        }
        else {
            //nodo malevolo
            probabilita_feedback_positivo = 0.5;
        }
    }

    double get_probabilità_feedback_positivo() {
        return this->probabilita_feedback_positivo;
    }
   
    void set_classe(int classe_disp) {
        classe = classe_disp;
    }

    int get_classe() {
        return this->classe;
    }

    void set_mu(int classe_disp) {
        if (classe_disp == 0) {
            //classe piu alta
            mu = 1.4;
        }
        else if (classe_disp == 1) {
            //classe media
            mu = 0.7;
        }
        else if (classe_disp == 2) {
            //classe peggiore
            mu = 0.025;
        }

    }


    double get_mu() {
        return this->mu;
    }
};

//costruttore
specifiche_nodo::specifiche_nodo() {
    this->id_nodo = 0;
    this->num_amici = 0;
    this->classe = 0;
    this->mu = 0;
    this->type_of_node = 0;
    this->probabilita_feedback_positivo = 0;
}

vector<specifiche_nodo> topologia;
vector<int> indici_amici_di_i;

tuple <int, int, int> stato;
map<tuple<int, int, int>, double> mapOfTuple;

int main() {
    //cosa voglio eseguire?
    int flag_prob_stato = 0;
    int flag_prob_perdita = 0;

    int indice_topologia = 0; //indice per il vettore di classi
    specifiche_nodo nodo_di_appoggio;
    specifiche_nodo reset;
    double variabile_appoggio = 0;
    string path = "SocialMatrix.txt";
    string path2 = "User_Info.txt";
    int contatore_amici = 0;
    string line;
    int social_index = 0;

    //SETTO LA TOPOLOGIA

    ifstream sim_file(path);
    ifstream sim_file_secondo(path2); //per benevolo malevolo e classe dispositivi
    for (indice_topologia = 0; indice_topologia < Number_of_nodes; indice_topologia++) {
        //leggere dato da file e mettere in appoggio
        
        nodo_di_appoggio.set_id_nodo(indice_topologia+1);
       
        if (sim_file.is_open()) {

            getline(sim_file, line);
            istringstream iss(line);
                for (social_index = 0; social_index < Number_of_nodes; social_index++)
                {
                    iss >> variabile_appoggio;
                   //cout << variabile_appoggio << '\n';
                    nodo_di_appoggio.S.push_back(variabile_appoggio);
                    if (variabile_appoggio > 0) 
                        contatore_amici++;
                }       
        }

        nodo_di_appoggio.set_num_amici(contatore_amici);
        contatore_amici = 0;

        if (sim_file_secondo.is_open()){

            getline(sim_file_secondo, );
            nodo_di_appoggio.set_type(0); //andrebbe preso anche questo dal simulatore
            nodo_di_appoggio.set_probabilità_feedback_positivo(nodo_di_appoggio.get_type());
            nodo_di_appoggio.set_classe(0); //andrebbe preso anche questo dal simulatore
            nodo_di_appoggio.set_mu(nodo_di_appoggio.get_classe());
        }
        topologia.push_back(nodo_di_appoggio);
        nodo_di_appoggio = reset;
    }
    sim_file.close();

    if (flag_prob_stato){
    
        stato = make_tuple(kj, Tj, nj);
        double P = 1;
        mapOfTuple[stato] = P;

        calcolo_denominatore_lambda_k();

        std::cout << "Stato: (" << kj << "," << Tj << "," << nj << ") Probabilita' di stato: " << P << endl;
        // int estract = mapOfTuple[make_tuple(i,j,k)];

        //la prima prob è uno, quindi inizio a calcolare dalla seconda
        nj++;

        for (Tj = 20; Tj <= max_passo; Tj++) {
            for (kj = 18; kj <= Tj - 2; kj++) {
                for (nj = 0; nj <= max_allocab_resources; nj++) {
                    if (Tj == 20 && kj == 18 && nj == 0)
                        continue;
                    std::cout << "Stato: (" << kj << "," << Tj << "," << nj << ") ";
                    stato = make_tuple(kj, Tj, nj);
                    P = state_probability(kj, Tj, nj);
                    mapOfTuple[stato] = P;
                    std::cout << "Probabilita' di stato: " << P << endl;

                }
            }
        }
    }
  
    if (flag_prob_perdita) {

        double prob_perdita = 0;
        prob_perdita = loss_probability();
    }

    return 0;
}

double state_probability(int kj, int Tj, int nj) {
   
    double probabilita = 0;
    double lambda_j = calcolo_lambda_j();
    
   
    //controllo sullo stato attraverso map

    double prob_prima_parte = 0;
    if (nj - 1 >= 0){
        prob_prima_parte = mapOfTuple[make_tuple(kj, Tj, nj - 1)] * lambda_j / (lambda_j + (nj - 1) * topologia[id_nodo_j-1].get_mu());
    }
    else{
        prob_prima_parte = 0;
    }
        
    double prob_seconda_parte = mapOfTuple[make_tuple(kj, Tj-1, nj+1)] * ((nj + 1) * topologia[id_nodo_j-1].get_mu() * (1-topologia[id_nodo_j - 1].get_probabilità_feedback_positivo())) / (lambda_j + (nj + 1) * topologia[id_nodo_j-1].get_mu());
    double prob_terza_parte = mapOfTuple[make_tuple(kj-1, Tj-1, nj+1)] * ((nj + 1) * topologia[id_nodo_j-1].get_mu() * topologia[id_nodo_j - 1].get_probabilità_feedback_positivo()) / (lambda_j + (nj + 1) * topologia[id_nodo_j-1].get_mu());

    probabilita = prob_prima_parte + prob_seconda_parte + prob_terza_parte;
    return probabilita;
}

void calcolo_denominatore_lambda_k() {
    int i = 0;
    for (i = 0; i < topologia.size(); i++) {
        denominatore_lambda_k = denominatore_lambda_k + topologia[i].get_probabilità_feedback_positivo() * topologia[i].get_num_amici();
    }
}
 
double calcolo_lambda_k(int id_nodo_k) {
    double lambda_k = 1;
    
    lambda_k = (lambda * topologia[id_nodo_k-1].get_probabilità_feedback_positivo() * topologia[id_nodo_k - 1].get_num_amici()) / denominatore_lambda_k;
    return lambda_k;
}

double calcolo_lambda_j() {
    double lambda_provider = 0;  //lambda_j risultato da ritornare
    vector<double> prob_richiesta_di_i_assegnata_a_j;
    
    int i;
  
    for (i = 0; i < Number_of_nodes; i++) {
        // dopo teorema delle prob totali sul numero di amici di i piu trusted di j
        if (topologia[id_nodo_j - 1].S[i] > 0) {  
            if (topologia[id_nodo_j - 1].S[i] * (double)(kj / Tj) > soglia) { //se sotto soglia inutile calcolarlo
                prob_richiesta_di_i_assegnata_a_j.push_back(trust_model(i));//passo indice id_amico di j
            }
            else
                prob_richiesta_di_i_assegnata_a_j.push_back(0);
        }    
        else
            prob_richiesta_di_i_assegnata_a_j.push_back(0);
        lambda_provider = lambda_provider + (lambda_ij * prob_richiesta_di_i_assegnata_a_j[i]);
    }
    return lambda_provider;
}

double trust_model(int indice_amico_di_j) {
    double valore_controllo_trust = 0;
    double probabilità_congiunta_k1 = 0;
    //double probabilità_congiunta_k2 = 0;
    double prodotto_prob_congiuta = 0;
    double prob_amico_piu_trusted_di_j = 1;
    double prob_amico_piu_trusted_di_j_parziale = 1;
    //double prob_amico_piu_trusted_di_j2 = 1;
    double prodotto_prob_amico_piu_trusted_di_j = 1;
    int amico_di_i = 0; // da verificare
    double lambda_k = 0;
    //double lambda_k2 = 0;
    int i;
    int j;
    int z;
    double T_km;
    int n_index=0;
    int n_0 = 0;
    double binomiale = 0;
   // int j2;
    specifiche_nodo km;

    for (i = 0; i < Number_of_nodes; i++)
    {
        if (topologia[indice_amico_di_j].S[i] > 0 && indice_amico_di_j != (id_nodo_j-1)) { 
               indici_amici_di_i.push_back(i);
        }
    }
            
    //for (i = 0; i < topologia[indice_amico_di_j].get_num_amici()-1; i++) {
    for (i = 0; i < 3; i++) {
        if (i == 0) {
           probabilità_congiunta_k1 = 1;
           prob_amico_piu_trusted_di_j = prob_binomiale(indici_amici_di_i, indice_amico_di_j); 
           valore_controllo_trust = probabilità_congiunta_k1 * prob_amico_piu_trusted_di_j;
        }
        else if(i == 1) {
            for (j = 0; j < indici_amici_di_i.size(); j++) {
                // equivale a dire che il k in questione non ha le risorse disponibili
                lambda_k = calcolo_lambda_k(indici_amici_di_i[j]);  //calcolo lambda k
                km = topologia[indici_amici_di_i[j]];
                //calcolo T_km
                if (km.get_type() == topologia[id_nodo_j - 1].get_type())
                    T_km = Tj;
                else {
                    T_km = proporzionalita_tk(km);
                }
                probabilità_congiunta_k1 = prob_di_blocco_generica(km, lambda_k); //questo sarà uguale alla prob di blocco di lambda k
                n_0 = (topologia[indice_amico_di_j].S[id_nodo_j - 1] * (double)(kj / Tj) * T_km) / km.S[indice_amico_di_j - 1];
                for (n_index = n_0 + 1; n_index <= T_km; n_index++) {
                    prob_amico_piu_trusted_di_j_parziale = prob_binomiale(indici_amici_di_i, indice_amico_di_j,j);//to check termini da passare perche devo passare tutto il vettore tranne indici_amici_di_i[j]
                    binomiale = factorial(T_km) / (factorial(n_index) * factorial(T_km - n_index));//calcolo binomiale tra T_km e n_index;
                    binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(), n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - n_index);
                    prob_amico_piu_trusted_di_j_parziale = prob_amico_piu_trusted_di_j_parziale * binomiale;
                    prob_amico_piu_trusted_di_j = prob_amico_piu_trusted_di_j + prob_amico_piu_trusted_di_j_parziale;
                }              
                valore_controllo_trust = valore_controllo_trust + (probabilità_congiunta_k1 * prob_amico_piu_trusted_di_j);
            }
        }
        /*
        caso uguale 2 solo se strettamente necessario
        else if (i == 2) {
            //to check here se posso fermarmi a due
            for (j = 0; j < amici_di_i.size(); j++) {
                lambda_k = calcolo_lambda_k(amici_di_i[j]);
                probabilità_congiunta_k1 = lambda_k; //aggiungere tutta la prob di blocco del kappesimo amico
                for (j2 = 0; j2 < amici_di_i.size(); j2++) {
                    if (amici_di_i[j] != amici_di_i[j2]) {
                        lambda_k2 = calcolo_lambda_k(amici_di_i[j2]);
                        probabilità_congiunta_k2 = lambda_k;
                        prodotto_prob_congiuta = probabilità_congiunta_k1 * probabilità_congiunta_k2;
                        //prob_amico_piu_trusted_di_j
                        //prob_amico_piu_trusted_di_j2
                        
                        //prodotto_prob_amico_piu_trusted_di_j = prob_amico_piu_trusted_di_j * prob_amico_piu_trusted_di_j2;
                        //valore_controllo_trust = valore_controllo_trust + (prodotto_prob_congiuta * prodotto_prob_amico_piu_trusted_di_j);
                    }
                }
            }               
        }
        */
    }
    return valore_controllo_trust;
}

double prob_binomiale(vector<int> indici_amici_di_i,int indice_amico_di_j) {

    double valore_probabilità_parziale = 1;
    double valore_binomiale_sommatoria_totale = 1;
    double binomiale = 0;
    int i = 0;
    int n_index = 0;
    int n0 = 0;
    int id_amico_di_i = 0;
    int T_km = 0;
    specifiche_nodo km;

    for (i = 0; i < indici_amici_di_i.size(); i++) {
        binomiale = 0;
        valore_binomiale_sommatoria_totale = 0;
        id_amico_di_i = indici_amici_di_i[i];
        km = topologia[id_amico_di_i];
        if (km.get_type() == topologia[id_nodo_j - 1].get_type())
            T_km = Tj;
        else {
            //to check: formula se non sono dello stesso tipo
            T_km = proporzionalita_tk(km);
        }
        //calcolare n0=(intero piu piccolo Sij*delta j dallo stato*T_km)/Sik;
        n0 = (topologia[indice_amico_di_j].S[id_nodo_j-1] * (double)(kj / Tj) * T_km) / km.S[indice_amico_di_j - 1];
        for (n_index = 0; n_index < n0; n_index++)
        {
            binomiale = factorial(T_km)/ (factorial(n_index) * factorial(T_km- n_index));//calcolo binomiale tra T_km e n_index;
            binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(),n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - n_index);
            valore_binomiale_sommatoria_totale = valore_binomiale_sommatoria_totale + binomiale;
        }
        valore_probabilità_parziale = valore_probabilità_parziale * valore_binomiale_sommatoria_totale;
    }

    return valore_probabilità_parziale;
}

double prob_binomiale(vector<int> indici_amici_di_i, int indice_amico_di_j, int indice_saltare) {

    double valore_probabilità_parziale = 1;
    double valore_binomiale_sommatoria_totale = 1;
    double binomiale = 0;
    int i = 0;
    int n_index = 0;
    int n0 = 0;
    int id_amico_di_i = 0;
    int T_km = 0;
    specifiche_nodo km;

    for (i = 0; i < indici_amici_di_i.size(); i++) {
        if(i != indice_saltare){
            binomiale = 0;
            valore_binomiale_sommatoria_totale = 0;
            id_amico_di_i = indici_amici_di_i[i];
            km = topologia[id_amico_di_i];
            if (km.get_type() == topologia[id_nodo_j - 1].get_type())
                T_km = Tj;
            else {
                //to check: formula se non sono dello stesso tipo
                T_km = proporzionalita_tk(km);
            }
            //calcolare n0=(intero piu piccolo Sij*delta j dallo stato*T_km)/Sik;
            n0 = (topologia[indice_amico_di_j].S[id_nodo_j - 1] * (double)(kj / Tj) * T_km) / km.S[indice_amico_di_j - 1];
            for (n_index = 0; n_index < n0; n_index++)
            {
                binomiale = factorial(T_km) / (factorial(n_index) * factorial(T_km - n_index));//calcolo binomiale tra T_km e n_index;
                binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(), n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - n_index);
                valore_binomiale_sommatoria_totale = valore_binomiale_sommatoria_totale + binomiale;
            }
            valore_probabilità_parziale = valore_probabilità_parziale * valore_binomiale_sommatoria_totale;
        }
    }

    return valore_probabilità_parziale;
}

int proporzionalita_tk(specifiche_nodo km) {
    int T_k = 0;

    //recuperare T_j dallo stato
    T_k = (Tj * km.get_probabilità_feedback_positivo() * km.get_num_amici()) / (topologia[id_nodo_j - 1].get_probabilità_feedback_positivo() * topologia[id_nodo_j - 1].get_num_amici());
    return T_k;

}

double prob_di_blocco_generica(specifiche_nodo km, double lambda_k) {
    double P_B_prima_parte = 1;
    double P_B_seconda_parte = 1;
    double P_B_terza_parte = 1;
    double sommatoria = 0;
    double risultato_sommatoria = 0;

    P_B_prima_parte = pow((lambda_k / km.get_mu()), (max_allocab_resources));
    P_B_seconda_parte = 1 / factorial(max_allocab_resources);


    int i;
    for (i = 1; i <= (max_allocab_resources); i++) {
        sommatoria = pow((lambda_k / km.get_mu()), i) * (1 / factorial(i));
        risultato_sommatoria = risultato_sommatoria + sommatoria;
    }
    P_B_terza_parte = 1 / risultato_sommatoria;

    double P_B = P_B_prima_parte * P_B_seconda_parte * P_B_terza_parte;

    return P_B;

}

//CASO NON OMOGENEO
/*double loss_probability() {
    double prob_perdita = 1;
    double traffico_perso = 0;
    int s = 0;

    cout << endl << "Probabilita' di perdita del sistema " << prob_perdita << endl;

    traffico_perso = lambda * prob_perdita;
    cout << endl << "traffico perso dal sistema " << traffico_perso << endl;
    return prob_perdita;
}*/



//CASO OMOGENEO
double loss_probability() {

    double prob_perdita = 1;
    double traffico_perso = 0;
    int s = 0;

    //assegno i valori dei tempi di servizio
    std::cout << endl << "Stampo valori di mu: " << endl;
    for (s = 0; s < num_classi_di_servizio; s++) {
        if (s == 0) {
            mu.push_back(1.4);
        }
        if (s == 1) {
            mu.push_back(0.7);
        }
        if (s >= 2) {
            mu.push_back(0.025);
        }
        std::cout << " mu[" << s+1 << "] = " << mu[s] << endl;
    }

    //assegno il numero di amici per ogni classe di servizio

    std::cout << endl << "Il numero totale di amici per object e': " << num_amici << endl << "Numero di amici per classe: " << endl;
    double proporz_classe1 = 9.3;
    double proporz_classe2 = 3.7;
    double proporz_classe3 = 0;
    double totale_proporzione = proporz_classe1 + proporz_classe2 + proporz_classe3;
    for (s = 0; s < num_classi_di_servizio; s++) {
        if (s == 0) {
            //num_amici_per_classe.push_back(round(num_amici * proporz_classe1 / totale_proporzione));
            num_amici_per_classe.push_back(num_amici * proporz_classe1 / totale_proporzione);
        }
        if (s == 1) {
            //num_amici_per_classe.push_back(round(num_amici * proporz_classe2 / totale_proporzione));
            num_amici_per_classe.push_back(num_amici * proporz_classe2 / totale_proporzione);
        }
        if (s == 2) {
           // num_amici_per_classe.push_back(round(num_amici * proporz_classe3 / totale_proporzione));
            num_amici_per_classe.push_back(num_amici * proporz_classe3 / totale_proporzione);
        }
       //num_amici_per_classe.push_back(num_amici/ num_classi_di_servizio);
        std::cout << " per classe[" << s+1 << "]: " << num_amici_per_classe[s] << endl;
    }

    Lambda.push_back(lambda / ceil(num_amici_per_classe[0]));
    //Lambda.push_back(lambda / 8);


    std::cout << endl << "Probabilita' di blocco per servente " << endl;
   // array_prob_blocco_per_classe.push_back(prob_blocco(0));
    //cout << " di classe[1]: " << array_prob_blocco_per_classe[0] << endl;
    for (s = 0; s < num_classi_di_servizio; s++) {
        array_prob_blocco_per_classe.push_back(prob_blocco(s));
        std::cout << " di classe[" << s+1 << "]: " << array_prob_blocco_per_classe[s] << endl;
        if (s < num_classi_di_servizio - 1)
            Lambda.push_back((lambda * array_prob_blocco_per_classe[s]) / ceil(num_amici_per_classe[s+1])); 
            //Lambda.push_back((lambda * array_prob_blocco_per_classe[s]) / 6);
    }
    
    
    for (s = 0; s < array_prob_blocco_per_classe.size(); s++) {
        prob_perdita = prob_perdita * array_prob_blocco_per_classe[s];
    }

    std::cout << endl << "Probabilita' di perdita del sistema "<< prob_perdita << endl;
   
    traffico_perso = lambda * prob_perdita;
    std::cout << endl << "traffico perso dal sistema " << traffico_perso << endl;
    return prob_perdita;
}

double prob_blocco(int index) {
    double P_B_prima_parte = 1;
    double P_B_seconda_parte = 1;
    double P_B_terza_parte = 1;
    double sommatoria = 0;
    double risultato_sommatoria = 0;

    P_B_prima_parte = pow((Lambda[index]/mu[index]),(max_allocab_resources));
    P_B_seconda_parte = 1/factorial(max_allocab_resources);


    int i;
    for (i = 1; i <= (max_allocab_resources); i++) {
        sommatoria = pow((Lambda[index] / mu[index]), i) * (1 / factorial(i));
        risultato_sommatoria = risultato_sommatoria + sommatoria;
    }
    P_B_terza_parte = 1/risultato_sommatoria;

    double P_B = P_B_prima_parte * P_B_seconda_parte * P_B_terza_parte;

    return P_B;
}

double factorial(int n) {
    if (n > 1)
        return n * factorial(n - 1);
    else
        return 1;
}

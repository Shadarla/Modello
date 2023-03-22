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


class specifiche_nodo
{
private:
    int id_nodo;
    int num_amici;
    int classe;
    double mu;
    bool type_of_node;
    double probabilita_feedback_positivo;
    int max_allocable_res;

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
            probabilita_feedback_positivo = 0.91;
        }
        else {
            //nodo malevolo
            probabilita_feedback_positivo = 0.77;
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
           mu = 1.42857;
            //mu = 1.2;
        }
        else if (classe_disp == 1) {
            //classe media
            mu = 0.71428;
            //mu = 0.6;
        }
        else if (classe_disp == 2) {
            //classe peggiore
            mu = 0.025;
        }

    }


    double get_mu() {
        return this->mu;
    }

    void set_maxallocableres(int classe_disp) {
        if (classe_disp == 0) {
            //classe piu alta
            max_allocable_res = 3;
            //max_allocable_res = 2;
        }
        else if (classe_disp == 1) {
            //classe media
            max_allocable_res = 2;
            //max_allocable_res = 1;
        }
        else if (classe_disp == 2) {
            //classe peggiore
            max_allocable_res = 1;
        }

    }

    int get_maxallocableres() {
        return this->max_allocable_res;
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
    this->max_allocable_res = 0;
}



//C:\Users\gianc\Documents\GitHub\SSIoT\Sim - n_services_1 - n_devices_25 - n_master_1 - lambda_10.000000 - tot_sim_500 - seed_3 - resource_ctrl_1 - qoe_ctrl_1

double state_probability ();
double calcolo_lambda_j();
double trust_model(int id_amico_di_j);
double calcolo_lambda_k (int id_nodo_k);
double loss_probability();
double prob_blocco(int index);
double factorial(int n);
double prob_binomiale(vector<int> amici_di_i, int id_amico_j);
double prob_binomiale(vector<int> amici_di_i, int id_amico_j, int j);
int proporzionalita_tk(specifiche_nodo km);
double prob_di_blocco_generica(specifiche_nodo km, double lambda_k);
double calcolo_max_passo();

int kj = 18;
int Tj = 20;
int nj = 0;
const int Number_of_nodes = 25;
const double lambda = 16;
//double lambda = 1;
const int num_classi_di_servizio = 2;
const double mu_j = 1.4; //verifica (potrebbe anche non servire piu)
const bool type_of_node = 0; //0 per benevolo, 1 per malevolo (potrebbe anche non servire piu)
const int max_allocab_resources = 2; //(potrebbe anche non servire piu)
const double alfa_soglia = 0.78;  //arbitrario
double max_passo = 1;
//int max_passo = Tj+9; //esiste una formuletta per calcolarlo
const int num_amici = 13; // supposto uguale per tutti (to check per caso non omogeneo e potrebbe non servire piu)
int id_nodo_j = 1; //id del nodo da valutare
const double lambda_i = (lambda / Number_of_nodes); //essendo omogenea la distribuzione al momento è uguale per tutti
const double soglia = 0.48;
//const double soglia = 0.08;


//vettori per prob blocco
vector<double> mu;
vector<double> num_amici_per_classe;
vector<double> array_prob_blocco_per_classe;
vector<double> Lambda;
double lambda_j = lambda_i; // alla prima iterazione imposto lambda_j uguale per tutti i provider
double lambda_k = lambda_i;
double denominatore_lambda_k = 1;



vector<specifiche_nodo> topologia;


tuple <int, int, int> stato;
map<tuple<int, int, int>, double> mapOfTuple;
map<tuple<int, int, int>, double> maplambda_j;

struct States {
    //vector<tuple <int, int, int>> state_param;
    map<tuple<int, int, int>, double> prob_map_of_tuple;
    map<tuple<int, int, int>, double> map_current_lambda_j;
    map<int, double> map_media;
};

vector<States> Markov_chains;

int main() {
    //cosa voglio eseguire?
    int flag_prob_stato = 1;
    int flag_prob_perdita = 0;
    int calcolomediaoutput = 1;
    

    int indice_topologia = 0; //indice per il vettore di classi
    specifiche_nodo nodo_di_appoggio;
    specifiche_nodo reset;
    double variabile_appoggio = 0;
    int variabile_appoggio_per_quel_fallito_di_antonio = 0;
    float appoggio_double=0;
    string path = "SocialMatrix.txt";
    string path2 = "1.UserInfo.txt";
    int contatore_amici = 0;
    string line;
    string line2;
    int social_index = 0;

    //********SETTO LA TOPOLOGIA*************//

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
            if (indice_topologia == 0) {
                getline(sim_file_secondo, line2);
            }
            getline(sim_file_secondo, line2);
            istringstream iss2(line2);
            for (int ind_ciclo = 0; ind_ciclo < 9; ind_ciclo++) {
                if (ind_ciclo == 6 || ind_ciclo == 7) {
                    iss2 >> appoggio_double;
                    //cout << appoggio_double << endl;
                }
                else {
                    iss2 >> variabile_appoggio_per_quel_fallito_di_antonio;
                    //cout << variabile_appoggio_per_quel_fallito_di_antonio << endl;
                    if (ind_ciclo == 4) {
                        //cout << variabile_appoggio_per_quel_fallito_di_antonio << endl;
                        if (variabile_appoggio_per_quel_fallito_di_antonio == 3)
                            variabile_appoggio_per_quel_fallito_di_antonio = 0;
                        else if (variabile_appoggio_per_quel_fallito_di_antonio == 2)
                            variabile_appoggio_per_quel_fallito_di_antonio = 1;
                        else if (variabile_appoggio_per_quel_fallito_di_antonio == 1)
                            variabile_appoggio_per_quel_fallito_di_antonio = 2;
                        //cout << variabile_appoggio_per_quel_fallito_di_antonio << endl;
                        nodo_di_appoggio.set_classe(variabile_appoggio_per_quel_fallito_di_antonio); 
                    }
                    if (ind_ciclo == 8) {
                        //cout << variabile_appoggio_per_quel_fallito_di_antonio << endl;
                        nodo_di_appoggio.set_type(variabile_appoggio_per_quel_fallito_di_antonio); 
                    }
                }
            }
            nodo_di_appoggio.set_probabilità_feedback_positivo(nodo_di_appoggio.get_type());
            nodo_di_appoggio.set_mu(nodo_di_appoggio.get_classe());
            nodo_di_appoggio.set_maxallocableres(nodo_di_appoggio.get_classe());
        }
        topologia.push_back(nodo_di_appoggio);
        nodo_di_appoggio = reset;
    }
    sim_file.close();
    sim_file_secondo.close();

    //************CALCOLO PROBABILITA' DI STATO**************//

    if (flag_prob_stato) {

            States reset_state;
            //max_passo = calcolo_max_passo();
            max_passo = 21;

            std::cout << "Calcolo max passo" << Tj + max_passo << endl;

            //for (id_nodo_j = 1; id_nodo_j <= 2; id_nodo_j++) {
            for (id_nodo_j = 1; id_nodo_j <= Number_of_nodes; id_nodo_j++) {


                States current_node_State = reset_state;

                kj = 18;
                Tj = 20;
                nj = 0;
                lambda_j = lambda_i;

                stato = make_tuple(kj, Tj, nj);
                //current_node_State.state_param.push_back(stato);
                double P = 1;
                mapOfTuple[stato] = P;
                current_node_State.prob_map_of_tuple[stato] = P;
                maplambda_j[stato] = lambda_j;
                current_node_State.map_current_lambda_j[stato] = lambda_j;
                denominatore_lambda_k = topologia[id_nodo_j - 1].get_probabilità_feedback_positivo() * (double)topologia[id_nodo_j - 1].get_num_amici();

                std::cout << "###### VALUTO IL NODO " << id_nodo_j << " di tipo: " << topologia[id_nodo_j - 1].get_probabilità_feedback_positivo() << "######" << endl;
                std::cout << "Stato: (" << kj << "," << Tj << "," << nj << ") Probabilita' di stato: " << P;
                std::cout << endl;
                // int estract = mapOfTuple[make_tuple(i,j,k)];

                //la prima prob è uno, quindi inizio a calcolare dalla seconda
                //nj++;

                for (Tj = 20; Tj <= 20 + max_passo; Tj++) {
                    for (kj = 18; kj <= Tj - 2; kj++) {
                        for (nj = 0; nj <= topologia[id_nodo_j - 1].get_maxallocableres(); nj++) { //controllare se minore o minore uguale
                            if (Tj == 20 && kj == 18 && nj == 0)
                                continue;
                            else {
                                std::cout << "Stato: (" << kj << "," << Tj << "," << nj << ") ";
                                stato = make_tuple(kj, Tj, nj);
                                P = state_probability();

                                // EVITARE NAN
                                if (!isnan(P)) {
                                    if (P == 0) {
                                        P = 0.00000000000001;
                                    }
                                }

                                else {
                                    P = 0.00000000000001;

                                }
                                mapOfTuple[stato] = P;

                                current_node_State.prob_map_of_tuple[stato] = P;
                                std::cout << "Probabilita' di stato: " << P;
                                std::cout << endl;
                                if (nj == 0) {
                                    lambda_j = calcolo_lambda_j();
                                    //EVITARE NAN
                                    if (!isnan(lambda_j)) {
                                        if (lambda_j == 0) {
                                            lambda_j = 0.00000000000001;
                                        }
                                    }
                                    else {
                                        lambda_j = 0.00000000000001;
                                    }
                                    maplambda_j[stato] = lambda_j;
                                }
                                else
                                {
                                    maplambda_j[stato] = maplambda_j[make_tuple(kj, Tj, nj - 1)];
                                }
                                current_node_State.map_current_lambda_j[stato] = maplambda_j[stato];
                                //std::cout << "lambda_j: " << lambda_j << endl;
                            }
                        }
                    }
                }
                std::cout << endl;
                std::cout << endl;

                double numeratore_media_output = 0;
                double sommatoria_media_output = 0;
                double media_output = 0;
                double check_prob_uguale_uno = 0;
                double check_prob_uguale_uno_appoggio = 0;
                double check_prob_app_prec = 0;


                if (calcolomediaoutput) {
                    //CALCOLO MEDIA DI OUTPUT
                    for (Tj = 20; Tj <= 20 + max_passo; Tj++) {
                        double weightdenom = 0;
                        for (kj = 18; kj <= Tj - 2; kj++) {
                            for (nj = 0; nj <= topologia[id_nodo_j - 1].get_maxallocableres(); nj++) {
                                weightdenom = weightdenom + mapOfTuple[make_tuple(kj, Tj, nj)];
                            }
                        }

                        for (kj = 18; kj <= Tj - 2; kj++) {
                            //per ogni k devo calcolare il peso 
                            media_output = (double)kj / (double)Tj;
                            for (nj = 0; nj <= topologia[id_nodo_j - 1].get_maxallocableres(); nj++) { //controllare se minore o minore uguale
                                //calcolare somma prob di stato
                                numeratore_media_output = numeratore_media_output + mapOfTuple[make_tuple(kj, Tj, nj)];


                            }

                            numeratore_media_output = (numeratore_media_output / weightdenom);
                            sommatoria_media_output = sommatoria_media_output + (media_output * numeratore_media_output);

                            numeratore_media_output = 0;
                        }


                        //ad ogni T devo salvare il fatto
                        std::cout << "Media Delta per Tj=: " << Tj << " : " << sommatoria_media_output << endl << endl;
                        //PRINT MEDIA OUTPUT
                        ofstream media_output;
                        media_output.open("mediaoutput.txt", ios::app);
                        if (media_output.is_open()) {
                            if (Tj == 20) {
                                media_output << "NODO " << id_nodo_j << ": \n";
                            }
                            media_output << "Media Delta per Tj=: \t" << Tj << "\t : \t" << sommatoria_media_output << "\n";
                            media_output.close();
                        }

                        current_node_State.map_media[Tj] = sommatoria_media_output;
                        sommatoria_media_output = 0;
                    }

                }
               

                // SALVARE L'INTERA CATENA IN UNA STRUTTURA DATI PER IL CALCOLO DEGLI OUTPUT
                Markov_chains.push_back(current_node_State);
            }

            
   


            // CALCOLO PROB DI BLOCCO DEL SISTEMA - VERSIONE AL MOMENTO FUNZIONANTE
            double prob_blocco_system = 0;
            double denominatore_lambda_totale_perdite = 0;


            //CALCOLO TRAFFICO PERSO VERSIONE DEL PROF

            int risorse_nodo = 0;

            for (id_nodo_j = 1; id_nodo_j <= Number_of_nodes; id_nodo_j++) {
                denominatore_lambda_totale_perdite = 0;
                Tj = 20 + ((int)(max_passo * ((double)topologia[id_nodo_j - 1].get_num_amici() / (double)topologia[0].get_num_amici()) * (topologia[id_nodo_j - 1].get_probabilità_feedback_positivo() / topologia[0].get_probabilità_feedback_positivo())));
                risorse_nodo = topologia[id_nodo_j - 1].get_maxallocableres();
                //ELIMINO NODI MALEVOLI?
                if (topologia[id_nodo_j - 1].get_type() == 0) {
                    for (kj = 18; kj <= Tj - 2; kj++) {
                        for (nj = 0; nj <= risorse_nodo; nj++) {
                            if (Markov_chains[id_nodo_j - 1].prob_map_of_tuple[make_tuple(kj, Tj, nj)] > 0) {
                                denominatore_lambda_totale_perdite = denominatore_lambda_totale_perdite + Markov_chains[id_nodo_j - 1].prob_map_of_tuple[make_tuple(kj, Tj, nj)];
                            }
                        }
                    }
                }
                if (topologia[id_nodo_j - 1].get_type() == 0) {
                    for (kj = 18; kj <= Tj - 2; kj++) {
                        //EVITARE NAN
                        if (((Markov_chains[id_nodo_j - 1].prob_map_of_tuple[make_tuple(kj, Tj, risorse_nodo)] * Markov_chains[id_nodo_j - 1].map_current_lambda_j[make_tuple(kj, Tj, risorse_nodo)]) / denominatore_lambda_totale_perdite) > 0)
                            prob_blocco_system = prob_blocco_system + ((Markov_chains[id_nodo_j - 1].prob_map_of_tuple[make_tuple(kj, Tj, risorse_nodo)] * Markov_chains[id_nodo_j - 1].map_current_lambda_j[make_tuple(kj, Tj, risorse_nodo)]) / denominatore_lambda_totale_perdite);
                    }
                }
            }


            cout << "Intensità di traffico perso usando i pesi come lambda: " << prob_blocco_system << endl;

            //PRINT INTENSITA TRAFFICO PERSO
            ofstream file_traffico;
            file_traffico.open("filetraffico.txt", ios::app);
            if (file_traffico.is_open()) {
                file_traffico << lambda << "\t" << prob_blocco_system << "\t" << 10 << "\n";
                file_traffico.close();
            }

            //CALCOLO PROBABILITA' CHE UN SERVICE PROVIDER DI CLASSE PIU ALTA ABBIA RISORSE LIBERE



            double high_class_av_num = 0;
            double high_class_av_den = 0;
            double prob_high_class_av = 0;

            for (id_nodo_j = 1; id_nodo_j <= Number_of_nodes; id_nodo_j++) {
                if (topologia[id_nodo_j - 1].get_maxallocableres() == 3) {
                    Tj = 20 + ((int)(max_passo * ((double)topologia[id_nodo_j - 1].get_num_amici() / (double)topologia[0].get_num_amici()) * (topologia[id_nodo_j - 1].get_probabilità_feedback_positivo() / topologia[0].get_probabilità_feedback_positivo())));

                    for (kj = 18; kj <= Tj; kj++) {
                        for (nj = 18; nj <= topologia[id_nodo_j-1].get_maxallocableres(); nj++) {
                            if (nj <= topologia[id_nodo_j - 1].get_maxallocableres() - 1) {
                                high_class_av_num = high_class_av_num + Markov_chains[id_nodo_j-1].prob_map_of_tuple[make_tuple(kj,Tj,nj)];
                            }
                            high_class_av_den = high_class_av_den + Markov_chains[id_nodo_j - 1].prob_map_of_tuple[make_tuple(kj, Tj, nj)];

                        }
                    }
                    prob_high_class_av = high_class_av_num/ high_class_av_den;
                    high_class_av_num = 0;
                    high_class_av_den = 0;

                    // STAMPA
                    std::cout << "Probabilità hc availability del nodo" << id_nodo_j << "per Tj=: " << Tj << " : " << prob_high_class_av << endl << endl;
                    ofstream highclass_output;
                    highclass_output.open("highclass_output.txt", ios::app);
                    if (highclass_output.is_open()) {
                        highclass_output << "NODO " << id_nodo_j << ": \n" << "Probabilità hc availability per Tj=: " << Tj << " : " << prob_high_class_av << ": \n" << ": \n";
                        highclass_output.close();
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

double state_probability() { 

    double probabilita = 0; //variabile da ritornare
    double prob_prima_parte = 0;
    double prob_seconda_parte = 0;
    double prob_terza_parte = 0;


    if (nj - 1 > 0) {
        // cout << endl << (lambda_j + (nj - 1)) << endl;
        //if (lambda_j == 0)
        //     prob_prima_parte = 0;
        // else
        //     prob_prima_parte = mapOfTuple[make_tuple(kj, Tj, nj - 1)] * lambda_j / (lambda_j + (nj - 1) * topologia[id_nodo_j-1].get_mu());
        prob_prima_parte = mapOfTuple[make_tuple(kj, Tj, nj - 1)] * maplambda_j[make_tuple(kj, Tj, nj - 1)] / (maplambda_j[make_tuple(kj, Tj, nj - 1)] + (double)(nj - 1) * topologia[id_nodo_j - 1].get_mu());
    }
    else if (nj - 1 == 0) {
        prob_prima_parte = mapOfTuple[make_tuple(kj, Tj, nj - 1)];
    }
    else {
        prob_prima_parte = 0;
    }

    if (nj == topologia[id_nodo_j - 1].get_maxallocableres()) {
        prob_seconda_parte = 0;
        prob_terza_parte = 0;
    }
    else if (nj + 1 == topologia[id_nodo_j - 1].get_maxallocableres()) {
        if (Tj - kj == 2) {
            prob_terza_parte = mapOfTuple[make_tuple(kj - 1, Tj - 1, nj + 1)] * topologia[id_nodo_j - 1].get_probabilità_feedback_positivo();
            prob_seconda_parte = 0;
        }
        else if (kj == 18) {
            prob_seconda_parte = mapOfTuple[make_tuple(kj, Tj - 1, nj + 1)] * (1 - topologia[id_nodo_j - 1].get_probabilità_feedback_positivo());
            prob_terza_parte = 0;
        }
        else {
            prob_seconda_parte = mapOfTuple[make_tuple(kj, Tj - 1, nj + 1)] * (1 - topologia[id_nodo_j - 1].get_probabilità_feedback_positivo());
            prob_terza_parte = mapOfTuple[make_tuple(kj - 1, Tj - 1, nj + 1)] * topologia[id_nodo_j - 1].get_probabilità_feedback_positivo();
        }

    }
    else {
        if (Tj - kj == 2) {
            prob_seconda_parte = 0;
            prob_terza_parte = mapOfTuple[make_tuple(kj - 1, Tj - 1, nj + 1)] * (double)(nj + 1) * (topologia[id_nodo_j - 1].get_mu()) * (topologia[id_nodo_j - 1].get_probabilità_feedback_positivo()) / (maplambda_j[make_tuple(kj - 1, Tj - 1, nj + 1)] + ((double)(nj + 1) * topologia[id_nodo_j - 1].get_mu()));
        }
        else if (kj == 18) {
            prob_seconda_parte = mapOfTuple[make_tuple(kj, Tj - 1, nj + 1)] * (double)(nj + 1) * topologia[id_nodo_j - 1].get_mu() * (1 - topologia[id_nodo_j - 1].get_probabilità_feedback_positivo()) / (maplambda_j[make_tuple(kj, Tj - 1, nj + 1)] + ((double)(nj + 1) * topologia[id_nodo_j - 1].get_mu()));
            prob_terza_parte = 0;
        }
        else {
            prob_seconda_parte = mapOfTuple[make_tuple(kj, Tj - 1, nj + 1)] * (double)(nj + 1) * topologia[id_nodo_j - 1].get_mu() * (1 - topologia[id_nodo_j - 1].get_probabilità_feedback_positivo()) / (maplambda_j[make_tuple(kj, Tj - 1, nj + 1)] + (double)(nj + 1) * topologia[id_nodo_j - 1].get_mu());
            prob_terza_parte = mapOfTuple[make_tuple(kj - 1, Tj - 1, nj + 1)] * (double)(nj + 1) * topologia[id_nodo_j - 1].get_mu() * topologia[id_nodo_j - 1].get_probabilità_feedback_positivo() / (maplambda_j[make_tuple(kj - 1, Tj - 1, nj + 1)] + (double)(nj + 1) * topologia[id_nodo_j - 1].get_mu());
        }
    
}

    

    //cout << endl << prob_prima_parte << " + " << prob_seconda_parte << " + " << prob_terza_parte << endl;

    probabilita = prob_prima_parte + prob_seconda_parte + prob_terza_parte;
    return probabilita;
}

double calcolo_lambda_j() {
    
    double lambda_provider = 0;  //lambda_j risultato da ritornare
    double prob_richiesta_di_i_assegnata_a_j = 0; 

    int i=0;

    for (i = 0; i < Number_of_nodes; i++) {
        // dopo teorema delle prob totali sul numero di amici di i piu trusted di j
        if (topologia[id_nodo_j - 1].S[i] > 0) {
        //if (topologia[id_nodo_j - 1].S[i] * ((double)kj / (double)Tj) >= soglia) { //se trust sotto soglia inutile calcolarlo SOGLIA!!
            // cout << (topologia[id_nodo_j - 1].S[i]) * ((double)kj / (double)Tj) << endl;
            prob_richiesta_di_i_assegnata_a_j = trust_model(i);//passo indice id_amico di j
            //prob_richiesta_di_i_assegnata_a_j = 0.5;//rigo di prova per debuggare
        }
        else{
            prob_richiesta_di_i_assegnata_a_j = 0;
        }
        lambda_provider = lambda_provider + (lambda_i * prob_richiesta_di_i_assegnata_a_j);
    }
    //cout << endl << lambda_provider << endl;
    return lambda_provider;
}

double trust_model(int indice_amico_di_j) { 
    //sto passando l'indice i di tutti i nodi che scorro nel for, in particolare quelli che sono amici soprasoglia di j
    double valore_controllo_trust = 0;
    double valore_controllo_trust_i_1 = 0;
    //double probabilità_congiunta_k1 = 0;
    double teta = 0;
    //double probabilità_congiunta_k2 = 0;
    double prodotto_prob_congiuta = 0;
    //double prob_amico_piu_trusted_di_j = 1;
    double omega = 1;
    double prob_amico_piu_trusted_di_j_parziale = 1;
    //double prob_amico_piu_trusted_di_j2 = 1;
    double prodotto_prob_amico_piu_trusted_di_j = 1;
    int amico_di_i = 0; // da verificare
    double lambda_k = lambda_i;
    //double lambda_k2 = 0;
    int i;
    //int j;
    int T_km=0 ;
    int n_index = 0;
    int n_0 = 0;
    double binomiale = 0;
    // int j2;
    specifiche_nodo km;
    vector<int> indici_amici_di_i;

    for (i = 0; i < Number_of_nodes; i++)
    {
        //per ogni amico di j che chiamo i, mi salvo a sua volta i suoi amici in un vettore tranne j
        //if (topologia[indice_amico_di_j].S[i] > 0 && indice_amico_di_j != (id_nodo_j - 1)) {
        if (topologia[indice_amico_di_j].S[i] > 0 && i != (id_nodo_j - 1)) {
            indici_amici_di_i.push_back(i);
        }
    }

    //CONFRONTO CHE PER I, I SUOI AMICI SIANO O MENO PIU TRUSTED DI J
    for (i = 0; i < topologia[indice_amico_di_j].get_num_amici()-1; i++) {
    //for (i = 0; i < 3; i++) {
        if (i == 0) {
            teta = 1;
            omega = 0;
            omega = prob_binomiale(indici_amici_di_i, indice_amico_di_j); //CONTROLLATA DOVREBBE STARE BENE
            valore_controllo_trust = teta * omega;
            if (!isnan(valore_controllo_trust)) {
                //cout << "Valore di controllo trust " << valore_controllo_trust << endl;
            }
            else {
                valore_controllo_trust = 0.00000000000001;
            }
            
        }
        else if (i == 1) {
            valore_controllo_trust_i_1 = 0;
            for (int j = 0; j < indici_amici_di_i.size(); j++) {
                //equivale a dire che il k in questione non ha le risorse disponibili
                lambda_k = calcolo_lambda_k(indici_amici_di_i[j]);  //calcolo lambda k-
                //lambda_k = lambda_i/5;
                km = topologia[indici_amici_di_i[j]];
                
                //calcolo T_km
                if (km.get_type() == topologia[id_nodo_j - 1].get_type())
                    T_km = Tj;
                else {
                    T_km = proporzionalita_tk(km);
                }

                teta = prob_di_blocco_generica(km, lambda_k); //questo sarà uguale alla prob di blocco di lambda k
                
                //EVITARE NAN
                if (teta > 0) {

                }
                else {
                    teta = 0.00000000000001;
                }
                
                omega = 0;
                n_0 = (int)floor((topologia[id_nodo_j - 1].S[indice_amico_di_j] * ((double)kj / (double)Tj) * T_km) / km.S[indice_amico_di_j]);
                //cout << km.S[indice_amico_di_j] << endl;
                //cout << n_0 << endl;
                
                for (n_index = n_0 + 1; n_index <= T_km; n_index++) {
                    //questa vorrebbe essere produttoria di sommatoria di elementi
                    prob_amico_piu_trusted_di_j_parziale = prob_binomiale(indici_amici_di_i, indice_amico_di_j, j);//to check termini da passare perche devo passare tutto il vettore tranne indici_amici_di_i[j]
                    //calcolo binomiale tra T_km e n_index;
                    binomiale = factorial(T_km) / (factorial(n_index) * factorial(T_km - n_index));
                    binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(), n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - (double)n_index);
                    prob_amico_piu_trusted_di_j_parziale = prob_amico_piu_trusted_di_j_parziale * binomiale;
                    
                    // EVITARE NAN
                    if (prob_amico_piu_trusted_di_j_parziale > 0) {
                        
                    }
                    else {
                        prob_amico_piu_trusted_di_j_parziale = 0.000000000001;
                    }
                    omega = omega + prob_amico_piu_trusted_di_j_parziale;

                    /*
                    //EVITARE NAN
                    if (omega > 0) {

                    }
                    else {
                        omega = 0.00000000000001;
                    }
                    */
                    
                }
                         

                valore_controllo_trust_i_1 = valore_controllo_trust_i_1 + (teta * omega);
                //valore_controllo_trust = valore_controllo_trust + (teta * omega * (double)indici_amici_di_i.size()); //QUESTA è UNA PROVA SOLO PER VEDERE SE PER PSI=2 CAMBIA DI MOLTO LA SITUA
            }
            valore_controllo_trust = valore_controllo_trust + valore_controllo_trust_i_1;
        }
        
        /*
        //QUESTA è UNA PROVA SOLO PER VEDERE SE PER PSI=2 CAMBIA DI MOLTO LA SITUA
        *         else
        {
            //EVITARE NAN
            if (valore_controllo_trust_i_1 > 0) {

            }
            else {
                valore_controllo_trust_i_1 = 0.00000000000001;
            }
            valore_controllo_trust = valore_controllo_trust + valore_controllo_trust_i_1;
        }
        */


    }

    //DEVO TOGLIERE IL NULL
    if (!isnan(valore_controllo_trust)) {
        //cout << "Valore di controllo trust " << valore_controllo_trust << endl;
    }
    else{
        valore_controllo_trust = 0.00000000000001;
    }
        
    return valore_controllo_trust;
}

double calcolo_max_passo() {
    //DEFINIRE LA FUNZIONE PER IL CALCOLO DEL MAX PASSO A SECONDA DEL RATE DI PERDITA CHE SI VUOLE PER LA SOGLIA
    double T_delta = 0;
    double denom_T_delta = 0;
    denom_T_delta = (alfa_soglia * (double)kj) - (0.5 * (double)Tj);
    //Tdelta= k0-alfa ko / (alfa ko/T0)-Pf ---- 0.5 è la probabilità di feedback positivo di un nodo malevolo al momento statica
    T_delta = ((double)kj - (alfa_soglia * (double)kj)) *(double)Tj / denom_T_delta;

    //return Tj + T_delta;
    return T_delta;
        
}

double calcolo_lambda_k(int indici_nodo_k) {
    lambda_k = lambda_i;
    double lambda_j_prec = 0;
    //qua va fatta la media degli stati da cui provengo con maplambda_j[maketuple()]
    //MEDIA LAMBDAJ PRECEDENTE
    
        if (kj == 18) {
            lambda_j_prec = maplambda_j[make_tuple(kj, Tj - 1, nj + 1)];
        }
        else if (Tj - kj == 2) {
            lambda_j_prec = maplambda_j[make_tuple(kj - 1, Tj - 1, nj + 1)];
        }
        else {
            lambda_j_prec = maplambda_j[make_tuple(kj, Tj - 1, nj + 1)] * ( mapOfTuple[make_tuple(kj, Tj - 1, nj + 1)] / ( mapOfTuple[make_tuple(kj, Tj - 1, nj + 1)] + mapOfTuple[make_tuple(kj - 1, Tj - 1, nj + 1)])) + maplambda_j[make_tuple(kj - 1, Tj - 1, nj + 1)] * ( mapOfTuple[make_tuple(kj - 1, Tj - 1, nj + 1)] / (mapOfTuple[make_tuple(kj, Tj - 1, nj + 1)] + mapOfTuple[make_tuple(kj - 1, Tj - 1, nj + 1)]));
        }
    
  

    //DEVO ELIMINARE IL NAN
    if (!isnan(lambda_j_prec)) {
        lambda_k = (lambda_j_prec * topologia[indici_nodo_k].get_probabilità_feedback_positivo() * topologia[indici_nodo_k].get_num_amici()) / denominatore_lambda_k;
    }
    else {
        lambda_j_prec = 0.00000000000001;
        lambda_k = (lambda_j_prec * topologia[indici_nodo_k].get_probabilità_feedback_positivo() * (double)topologia[indici_nodo_k].get_num_amici()) / denominatore_lambda_k;
    }
    //std::cout << endl;
    //std::cout << "lambda_k del nodo " << topologia[indici_nodo_k].get_id_nodo() << ": " << lambda_k << endl;
    return lambda_k;
}

double prob_binomiale(vector<int> indici_amici_di_i,int indice_amico_di_j) {

    double valore_probabilità_parziale = 1;
    double valore_binomiale_sommatoria_totale = 0;
    double binomiale = 0;
    int i = 0;
    int n_index = 0;
    int n0 = 0;
    int id_amico_di_i = 0;
    int T_km = 0;
    specifiche_nodo km;

    for (i = 0; i < indici_amici_di_i.size(); i++) {
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
        n0 = (int)floor((topologia[id_nodo_j-1].S[indice_amico_di_j] * ((double)kj / (double)Tj) * (double)T_km) / km.S[indice_amico_di_j]);

        for (n_index = 0; n_index < n0; n_index++)
        {
            binomiale = factorial(T_km)/ (factorial(n_index) * factorial(T_km- n_index));//calcolo binomiale tra T_km e n_index;
            binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(),n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - n_index);
            valore_binomiale_sommatoria_totale = valore_binomiale_sommatoria_totale + binomiale;
        }
        valore_probabilità_parziale = valore_probabilità_parziale * valore_binomiale_sommatoria_totale;
    }
    //DEVO ELIMINARE IL NAN
    if (!isnan(valore_probabilità_parziale)) {
       
    }
    else {
        valore_probabilità_parziale = 0.00000000000001;
    }

    return valore_probabilità_parziale;
}

double prob_binomiale(vector<int> indici_amici_di_i, int indice_amico_di_j, int indice_saltare) {

    double valore_probabilità_parziale = 1;
    double valore_binomiale_sommatoria_totale = 0;
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
                //formula se non sono dello stesso tipo
                T_km = proporzionalita_tk(km);
            }

            //calcolare n0=(intero piu piccolo Sij*delta j dallo stato*T_km)/Sik;
            n0 = (int)floor((topologia[indice_amico_di_j].S[id_nodo_j - 1] * ((double)kj / (double)Tj) * T_km) / km.S[indice_amico_di_j]);
            for (n_index = 0; n_index < n0; n_index++)
            {
                binomiale = factorial(T_km) / (factorial(n_index) * factorial(T_km - n_index));//calcolo binomiale tra T_km e n_index;
                binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(), n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - n_index);
                //DEVO ELIMINARE IL NAN
                if (!isnan(binomiale)) {

                }
                else {
                    binomiale = 0.00000000000001;
                }
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
    T_k = ((int)((double)Tj * km.get_probabilità_feedback_positivo() * (double)km.get_num_amici() / (topologia[id_nodo_j - 1].get_probabilità_feedback_positivo() * (double)topologia[id_nodo_j - 1].get_num_amici())));
    return T_k;

}

double prob_di_blocco_generica(specifiche_nodo km, double lambda_k) {
    double P_B_prima_parte = 1;
    double P_B_seconda_parte = 1;
    double P_B_terza_parte = 1;
    double sommatoria = 0;
    double risultato_sommatoria = 0;

    P_B_prima_parte = pow((lambda_k / km.get_mu()), (double)(km.get_maxallocableres()));
    P_B_seconda_parte = 1 / factorial(km.get_maxallocableres());


    int i;
    for (i = 1; i <= (km.get_maxallocableres()); i++) {
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

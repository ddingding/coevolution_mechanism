bayes_hierarchical_syn = """
data{
    //replicate 1
    // synonymous count data
    int<lower=0> K;             // number of synonymous mutants
    int<lower=0> c_pre1[K];     // counts synonymous before selection
    int<lower=0> c_aft1[K];     // counts synonymous after selection

    //replicate 2
    // synonymous count data
    int<lower=0> c_pre2[K];     // counts synonymous before selection
    int<lower=0> c_aft2[K];     // counts synonymous after selection

    // mutant data.
    int<lower=0> Km;            // number of observed aa mutants
    int<lower=0> N;             // number of total codon observations
    int s[N];                   //group number of codon belonging to which amino acid group
                                // [1 1 2 2 2 3 3 3 3 ...] means that first 2 codons belong to amino acid group 1,
                                // next 3 codons  belong to amino acid group 2

    // replicate 1 mutants
    int<lower=0> c_pre_m1[N];   //all codon observations pre
    int<lower=0> c_aft_m1[N];

    // rep 2 mutants
    int<lower=0> c_pre_m2[N];   //all codon observations pre
    int<lower=0> c_aft_m2[N];
}

parameters {
    // growth rates w live on the reals
    // growth ratios like exp(w) live on the postive reals, and are modelled as a lognormal distribution
    //

    //global parameters between both synonymous and mutant codon observations
    vector[K] b1;                   // log rate parameter of c_pre for replicate 1,
                                    // shared between synonymous and mutant observations
    vector[K] b2;                   // log rate parameter of c_pre for replicate 2
    real<lower=0> sigma_b1;         // variance for lognormal model for rate of counts pre-selection, replicate 1
    real<lower=0> sigma_b2;         // analogous for replicate 2
    real rep_scale;                 // scale parameter for replicate 1 to replicate 2 differences in
                                    // observed growth rate differences. This is due to differences
                                    // in sequencing depth between replicates
    real<lower=0> sigma_rep;        // shared variance in codon level growth rates between replicates

    real<lower=0> sigma_w;          // variance of growth rates w and w_cod_m drawn from N(0,sigma_w)



    //synonymous params
    real w_aa;                         // the synonymous group growth rate (lives in log space of ratio),
                                    // or mean of lognormal for growth ratios
    vector[K] w_cod_std;            //  reparameterization of w_cod[k] ~ N(w, sigma_w_cod_std),
                                    // the growth rates for each individual synonymous mutant observation
    real<lower=0> sigma_w_cod_std;  // shared variance for w_cod
    vector[K] w_cod_r1_std;         // reparameterized to get replicate 1 growth rate, w_cod_r1 ~ N(w_cod, sigma_rep)
    vector[K] w_cod_r2_std;         // stdev for rep 2


    // mutants params
    vector[N] b_m1;                         // mutant log rate parameter for replicate 1 c_pre1
    vector[N] b_m2;                         // mutant log rate parameter for replicate 2 c_pre2
    vector[Km] w_aa_m;                         // amino acid level growth rate
    vector[N] w_cod_m_std;                  // reparameterized normal growth rate for each codon,
                                            // where synonymous codons are drawn from one normal w_cod_m[k]
                                            // w_cod_m[n] ~ N(w_cod_m[k], sigma_w_cod_m_std)
    vector<lower=0>[Km] sigma_w_cod_m_std;  // variance per amino acid group
    vector<lower=0>[Km] sigma_rep_m;
    vector[N] w_cod_m_r1_std;               // reparametrized normal growth rate for codon n and replicate 1:
                                            // w_cod_m_r1[n] ~ N(w_cod_m[n], sigma_rep)
    vector[N] w_cod_m_r2_std;               // as above for replicate 2


}

transformed parameters {
    // for synonymous
    vector<lower=0>[K] f_pre1;              // rate parameter poisson emitting counts pre-selection replicate 1
    vector<lower=0>[K] f_pre2;              // cf above for rep2
    vector<lower=0>[K] f_aft1;              // rate parameter for counts after selection replicate 1
    vector<lower=0>[K] f_aft2;              // cf above for rep2

    vector[K] w_cod;                        // parameter for each synonymous mutant
    vector[K] w_cod_r1;                     // synonymous growth rates for each K observation, replicate 1
    vector[K] w_cod_r2;                     // for rep2


    // for mutants
    vector<lower=0>[N] f_pre_m1;            // for the single mutants
    vector<lower=0>[N] f_pre_m2;
    vector<lower=0>[N] f_aft_m1;
    vector<lower=0>[N] f_aft_m2;

    vector[N] w_cod_m;                      // individual mutant growth rates for each codon
    vector[N] w_cod_m_r1;                   // individual mutant growht rates for each codon, replicate 1
    vector[N] w_cod_m_r2;


    for (k in 1:K){
        w_cod[k] = w_aa + w_cod_std[k] * sigma_w_cod_std;           // reparameterized w_cod[k] ~ Normal(w, sigma_w_cod)
        w_cod_r1[k] = w_cod[k] + w_cod_r1_std[k] * sigma_rep/10;       // reparameterized:
                                                                    // w_cod_r1[k] ~ Normal(w_cod[k], sigma_rep)
        w_cod_r2[k] = w_cod[k] + w_cod_r2_std[k] * sigma_rep/10;       // analogous for rep 2
    }

    //synonymous
    f_pre1 = exp(b1);                           // transforming b1 from reals to postive reals,
                                                // implying  f_pre1 ~ lognormal(0, sigma_b1)
    f_aft1 = exp(b1 + w_cod_r1);                // transforming w_cod_r1 from reals to positive reals,
                                                // implying f_aft1 = exp(b1) * r
                                                // r ~ lognormal(w_cod, sigma_rep)
    f_pre2 = exp(b2);                           // replicate 2
    f_aft2 = exp(b2 + w_cod_r2 + rep_scale);    // scaling the replicate 2 growth rates by one global scaling factor


    //per mutant
    for (n in 1:N){
        w_cod_m[n] = w_aa_m[s[n]] + w_cod_m_std[n] * sigma_w_cod_m_std[s[n]];   // reparameterized
                                                                                // w_cod_m[n] ~ N(w_aa_m[k],
                                                                                // sigma_w_cod_m[k])
        // draw 2 independent codon observations
        w_cod_m_r1[n] = w_cod_m[n] + w_cod_m_r1_std[n] *  sigma_rep_m[s[n]]/10;            // reparameterized:
                                                                                // w_cod_m_r1[n] ~ N(w_cod_m[n],
                                                                                // sigma_rep)
        w_cod_m_r2[n] = w_cod_m[n] + w_cod_m_r2_std[n] *  sigma_rep_m[s[n]]/10;            // analogous for rep2

        f_pre_m1[n] = exp(b_m1[n]);                             // implies f_pre_m1[n] ~ lognormal(0, sigma_b1)
        f_aft_m1[n] = exp(b_m1[n] + w_cod_m_r1[n]);             // transforming w_cod_r1 from reals to positive reals
                                                                // implying f_aft1 = exp(b1) * r
        f_pre_m2[n] = exp(b_m2[n]);
        f_aft_m2[n] = exp(b_m2[n] + w_cod_m_r2[n] + rep_scale);
    }
}


model {
    // correlation between replicates
    sigma_rep ~ exponential(3);                 // prior belief that replicate sigma is low
    //sigma_rep ~ cauchy(0,0.3);

    // for the synonymous counts
    w_aa ~ normal(-2,5);                        // just sample this from a pretty vague prior
    w_cod_std ~ normal(0,1);                    // implies: w_cod[k] ~ normal(w, sigma_w_cod_std)
    sigma_w_cod_std ~ exponential(8);           // prior belief that w_cods are drawn from a tight distribution
    w_cod_r1_std ~ normal(0,1);                 // implies: w_cod_r1[k] ~ normal(w_cod[k], sigma_rep;)
    w_cod_r2_std ~ normal(0,1);                 // implies: w_cod_r2[k] ~ normal(w_cod[k], sigma_rep;)

    // synonymous  muts
    for (i in 1:K) {
        b1[i] ~ normal(0,sigma_b1);             // draw the log rate parameter for f_pre1
        c_pre1[i] ~ poisson(f_pre1[i]);         // emitting poissons for pre and after selection counts
        c_aft1[i] ~ poisson(f_aft1[i]);

        b2[i] ~ normal(0,sigma_b2);             // cf rep 2
        c_pre2[i] ~ poisson(f_pre2[i]);
        c_aft2[i] ~ poisson(f_aft2[i]);
    }

    // amino acid level group growth rates
    for (k in 1:Km){
        w_aa_m[k] ~ normal(-2,2);             // draw the amino acid level group growth rates
        sigma_w_cod_m_std[k] ~ exponential(8);  // prior belief that this is small
                                                // each amino acid group has its own stdev
        sigma_rep_m[k] ~ exponential(3);
    }

    // codon level growth ratea
    for (n in 1:N){
        // for each codon mutation sample a random variable
        w_cod_m_std[n] ~ normal(0,1);           // implies w_cod_m[n] ~ N(w_aa_m[k], sigma_w_aa_m[k])
        //for individual replicate muts
        w_cod_m_r1_std[n] ~ normal(0,1);            // implies: w_cod_m_r1[k] ~ normal(w_cod_m[n], sigma_rep)
        w_cod_m_r2_std[n] ~ normal(0,1);            // implies: w_cod_m_r2[k] ~ normal(w_cod_m[n], sigma_rep)
        //
        b_m1[n] ~ normal(0,sigma_b1);           // drawing log rate parameters as before
        b_m2[n] ~ normal(0,sigma_b2);
        c_pre_m1[n] ~ poisson(f_pre_m1[n]);
        c_pre_m2[n] ~ poisson(f_pre_m2[n]);
        c_aft_m1[n] ~ poisson(f_aft_m1[n]);
        c_aft_m2[n] ~ poisson(f_aft_m2[n]);
    }


}

generated quantities{

    // quantities of interest
    vector[Km] diff_w_aa_m; // difference in growth rates, which live on reals
    vector[Km] diff_r; // difference in growth ratios, which live on positive reals
    vector[Km] ratio_r; // ratio of growth ratio

    // for PPC
    // replicate 1
    //synonymous counts
    int<lower=0> c_pre1_pred[K];
    int<lower=0> c_aft1_pred[K];
    //mutant counts
    int<lower=0> c_pre_m1_pred[N];
    int<lower=0> c_aft_m1_pred[N];

    // replicate 2
    //synonymous counts
    int<lower=0> c_pre2_pred[K];
    int<lower=0> c_aft2_pred[K];
    //mutant counts
    int<lower=0> c_pre_m2_pred[N];
    int<lower=0> c_aft_m2_pred[N];

    // quantities of interest
    diff_r = exp(w_aa_m) - exp(w_aa);
    ratio_r = exp(w_aa_m) / exp(w_aa);
    diff_w_aa_m = w_aa_m - w_aa;

    // generate synonymous pre and post counts
    for (i in 1:K) {
        //rep1
        c_pre1_pred[i] = poisson_rng(f_pre1[i]);
        c_aft1_pred[i] = poisson_rng(f_aft1[i]);
        // rep2
        c_pre2_pred[i] = poisson_rng(f_pre2[i]);
        c_aft2_pred[i] = poisson_rng(f_aft2[i]);
    }

    for (i in 1:N) {
        //rep1
        c_pre_m1_pred[i] = poisson_rng(f_pre_m1[i]);
        c_aft_m1_pred[i] = poisson_rng(f_aft_m1[i]);
        //rep2
        c_pre_m2_pred[i] = poisson_rng(f_pre_m2[i]);
        c_aft_m2_pred[i] = poisson_rng(f_aft_m2[i]);
    }

}

"""

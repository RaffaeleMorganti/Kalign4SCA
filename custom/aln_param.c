#include "tldevel.h"

#include "kalign/kalign.h"
#include "msa_struct.h"

#define ALN_PARAM_IMPORT
#include "aln_param.h"

static int set_subm_gaps_protein(struct aln_param *ap);
static int set_subm_gaps_DNA(struct aln_param *ap);

int aln_param_init(struct aln_param **aln_param,int biotype , int n_threads, int type, float gpo, float gpe, float tgpe)
{
        struct aln_param* ap = NULL;

        /* Allocate  */
        MMALLOC(ap, sizeof(struct aln_param));
        ap->subm = NULL;

        ap->nthreads = n_threads;
        MMALLOC(ap->subm,sizeof (float*) * 23);

        for (int i = 23;i--;){
                ap->subm[i] = NULL;
                MMALLOC(ap->subm[i],sizeof(float) * 23);
                for (int j = 23;j--;){
                        ap->subm[i][j] = 0.0f;
                }
        }
        if(biotype == ALN_BIOTYPE_DNA){
                /* include/kalign/ */
                switch (type) {
                case KALIGN_TYPE_PROTEIN:
                        ERROR_MSG("Detected DNA sequences but --type protein option was selected.");
                        break;
                case KALIGN_TYPE_DNA:
                case KALIGN_TYPE_DNA_INTERNAL:
                case KALIGN_TYPE_RNA:
                default:
                        set_subm_gaps_DNA(ap);
                        break;
                }
        }else if(biotype == ALN_BIOTYPE_PROTEIN){
                switch (type) {
                case KALIGN_TYPE_DNA:
                        ERROR_MSG("Detected protein sequences but --type dna option was selected.");
                        break;
                case KALIGN_TYPE_DNA_INTERNAL:
                        ERROR_MSG("Detected protein sequences but --type internal  option was selected.");
                        break;
                case KALIGN_TYPE_RNA:
                        ERROR_MSG("Detected protein sequences but --type rna option was selected.");
                        break;
                case KALIGN_TYPE_PROTEIN:
                case KALIGN_TYPE_PROTEIN_DIVERGENT:
                default:
                        set_subm_gaps_protein(ap);
                        /* set_subm_gaps_gon250(ap); */
                        break;
                }
        }else{
                ERROR_MSG("Unable to determine what alphabet to use.");
        }

        if(gpo >= 0.0){
                ap->gpo = gpo;
        }
        if(gpe >= 0.0){
                ap->gpe = gpe;
        }

        if(gpe >= 0.0){
                ap->tgpe = tgpe;
        }
        /* LOG_MSG("%f %f %f", ap->gpo, ap->gpe, ap->tgpe); */
        *aln_param = ap;
        return OK;
ERROR:
        aln_param_free(ap);
        return FAIL;
}



int set_subm_gaps_protein(struct aln_param* ap)
{
        // ABC DEF GHI JKL MNO PQR
        float score[18] = {31,30,28,24,16,0,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29,-30,-31};
        for (int i = 0; i < 18; i++) {
                for (int j = 0; j <= i; j++) {
                        ap->subm[i][j] = score[i-j];
                        ap->subm[j][i] = score[i-j];
                }
        }

        ap->gpo = 20;
        ap->gpe = 10;
        ap->tgpe = 1;
        return OK;
}

int set_subm_gaps_DNA(struct aln_param *ap)
{
        // UVWXYZ
        float score[6] = {4,2,0,0,0,2};
        for (int i = 0; i < 6; i++) {
                for (int j = 0; j <= i; j++) {
                        ap->subm[i][j] = score[i-j];
                        ap->subm[j][i] = score[i-j];
                }
        }

        ap->gpo = 2;
        ap->gpe = 1;
        ap->tgpe = 0;
        return OK;
}


void aln_param_free(struct aln_param *ap)
{
        if(ap){
                if(ap->subm){
                        for (int i = 23;i--;){
                                MFREE(ap->subm[i]);
                        }
                        MFREE(ap->subm);
                }
                MFREE(ap);
        }
}

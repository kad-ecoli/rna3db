/* parameters for RNA base pairing */
#ifndef BPstat_HPP
#define BPstat_HPP 1
const float Pm_mu    = 70.036; const float Pm_sd    =53.306; // P[i-1]-P[i]-P[j]-P[j-1]
const float Pp_mu    = 81.983; const float Pp_sd    =24.401; // P[i+1]-P[i]-P[j]-P[j+1]
const float O5m_mu   = 75.236; const float O5m_sd   =52.434; // O5'[i-1]-O5'[i]-O5'[j]-O5'[j-1]
const float O5p_mu   = 96.081; const float O5p_sd   =21.428; // O5'[i+1]-O5'[i]-O5'[j]-O5'[j+1]
const float C5m_mu   = 55.274; const float C5m_sd   =51.499; // C5'[i-1]-C5'[i]-C5'[j]-C5'[j-1]
const float C5p_mu   = 99.530; const float C5p_sd   =24.674; // C5'[i+1]-C5'[i]-C5'[j]-C5'[j+1]
const float C4m_mu   = 61.550; const float C4m_sd   =48.747; // C4'[i-1]-C4'[i]-C4'[j]-C4'[j-1]
const float C4p_mu   =114.858; const float C4p_sd   =25.515; // C4'[i+1]-C4'[i]-C4'[j]-C4'[j+1]
const float C3m_mu   = 83.423; const float C3m_sd   =31.092; // C3'[i-1]-C3'[i]-C3'[j]-C3'[j-1]
const float C3p_mu   =132.495; const float C3p_sd   =28.286; // C3'[i+1]-C3'[i]-C3'[j]-C3'[j+1]
const float C2m_mu   = 92.312; const float C2m_sd   =35.009; // C2'[i-1]-C2'[i]-C2'[j]-C2'[j-1]
const float C2p_mu   =146.662; const float C2p_sd   =31.430; // C2'[i+1]-C2'[i]-C2'[j]-C2'[j+1]
const float C1m_mu   = 70.909; const float C1m_sd   =54.285; // C1'[i-1]-C1'[i]-C1'[j]-C1'[j-1]
const float C1p_mu   =133.920; const float C1p_sd   =29.379; // C1'[i+1]-C1'[i]-C1'[j]-C1'[j+1]
const float O4m_mu   = 59.676; const float O4m_sd   =53.181; // O4'[i-1]-O4'[i]-O4'[j]-O4'[j-1]
const float O4p_mu   =113.180; const float O4p_sd   =25.776; // O4'[i+1]-O4'[i]-O4'[j]-O4'[j+1]
const float O3m_mu   = 82.379; const float O3m_sd   =26.246; // O3'[i-1]-O3'[i]-O3'[j]-O3'[j-1]
const float O3p_mu   =136.875; const float O3p_sd   =30.475; // O3'[i+1]-O3'[i]-O3'[j]-O3'[j+1]
const float P44P_mu  = -0.990; const float P44P_sd  =33.603; // P[i]-C4'[i]-C4'[j]-P[j]
const float C4114C_mu=172.313; const float C4114C_sd=43.508; // C4'[i]-C1'[i]-C1'[j]-C4'[j]
const float PP_mu    = 18.519; const float PP_sd    = 1.295; // P[i]-P[j]
const float O5O5_mu  = 16.806; const float O5O5_sd  = 0.790; // O5'[i]-O5'[j]
const float C5C5_mu  = 17.160; const float C5C5_sd  = 0.518; // C5'[i]-C5'[j]
const float C4C4_mu  = 15.018; const float C4C4_sd  = 0.419; // C4'[i]-C4'[j]
const float C3C3_mu  = 13.676; const float C3C3_sd  = 0.498; // C3'[i]-C3'[j]
const float C2C2_mu  = 11.046; const float C2C2_sd  = 0.515; // C2'[i]-C2'[j]
const float C1C1_mu  = 10.659; const float C1C1_sd  = 0.336; // C1'[i]-C1'[j]
const float O4O4_mu  = 13.232; const float O4O4_sd  = 0.399; // O4'[i]-O4'[j]
const float O3O3_mu  = 15.265; const float O3O3_sd  = 0.645; // O3'[i]-O3'[j]
const float NN_mu    =  8.928; const float NN_sd    = 0.253; // N[i]-N[j]
const float aPm_mu   =104.731; const float aPm_sd   =50.284; // <P[i-1]P[i],P[j+1]P[j]>
const float aPp_mu   =107.925; const float aPp_sd   =32.776; // <P[i+1]P[i],P[j-1]P[j]>
const float aO5m_mu  = 94.387; const float aO5m_sd  =42.182; // <O5'[i-1]O5'[i],O5'[j+1]O5'[j]>
const float aO5p_mu  = 97.722; const float aO5p_sd  =19.370; // <O5'[i+1]O5'[i],O5'[j-1]O5'[j]>
const float aC5m_mu  = 91.152; const float aC5m_sd  =42.549; // <C5'[i-1]C5'[i],C5'[j+1]C5'[j]>
const float aC5p_mu  = 94.525; const float aC5p_sd  =21.115; // <C5'[i+1]C5'[i],C5'[j-1]C5'[j]>
const float aC4m_mu  = 78.536; const float aC4m_sd  =40.484; // <C4'[i-1]C4'[i],C4'[j+1]C4'[j]>
const float aC4p_mu  = 81.262; const float aC4p_sd  =19.154; // <C4'[i+1]C4'[i],C4'[j-1]C4'[j]>
const float aC3m_mu  = 68.348; const float aC3m_sd  =39.923; // <C3'[i-1]C3'[i],C3'[j+1]C3'[j]>
const float aC3p_mu  = 70.502; const float aC3p_sd  =19.120; // <C3'[i+1]C3'[i],C3'[j-1]C3'[j]>
const float aC2m_mu  = 56.490; const float aC2m_sd  =40.294; // <C2'[i-1]C2'[i],C2'[j+1]C2'[j]>
const float aC2p_mu  = 58.240; const float aC2p_sd  =21.273; // <C2'[i+1]C2'[i],C2'[j-1]C2'[j]>
const float aC1m_mu  = 61.593; const float aC1m_sd  =40.418; // <C1'[i-1]C1'[i],C1'[j+1]C1'[j]>
const float aC1p_mu  = 63.502; const float aC1p_sd  =21.469; // <C1'[i+1]C1'[i],C1'[j-1]C1'[j]>
const float aO4m_mu  = 75.154; const float aO4m_sd  =40.504; // <O4'[i-1]O4'[i],O4'[j+1]O4'[j]>
const float aO4p_mu  = 77.562; const float aO4p_sd  =20.099; // <O4'[i+1]O4'[i],O4'[j-1]O4'[j]>
const float aO3m_mu  = 67.558; const float aO3m_sd  =39.669; // <O3'[i-1]O3'[i],O3'[j+1]O3'[j]>
const float aO3p_mu  = 69.837; const float aO3p_sd  =18.743; // <O3'[i+1]O3'[i],O3'[j-1]O3'[j]>
const float aPC_mu   = 58.803; const float aPC_sd   =27.207; // <P[i]C4'[i],P[j]C4'[j]>
const float aCC_mu   =165.213; const float aCC_sd   =13.607; // <C4'[i]C1'[i],C4'[j]C1'[j]>
#endif

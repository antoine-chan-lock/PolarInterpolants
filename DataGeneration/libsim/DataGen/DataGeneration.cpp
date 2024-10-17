#include "DataGeneration.h"
#include <igl/find.h>
#include "progressbar.hpp"

DataGeneration::DataGeneration()
{
}

DataGeneration::~DataGeneration()
{
}

DataGeneration::DataGeneration(const string &matid)
{
    m_matid = matid;
    m_matdir = "./data/microstructures/mat"+matid;
    m_slv = Solver(m_matdir);
}

void DataGeneration::sampling()
{
    int density = 8;
    ArrayXd tmp = ArrayXd::LinSpaced(density+1, Tools::deg2rad(0),Tools::deg2rad(180));

    ArrayXd m_va = tmp(seq(0,tmp.size()-2));
    ArrayXd m_vs = ArrayXd::LinSpaced(16, -0.5, 1);
    ArrayXd m_vktmp = ArrayXd::LinSpaced(16, 0, 1);
    ArrayXd m_vk(m_vktmp.size()*2); 
    m_vk << m_vktmp, -m_vktmp;
    ArrayXd m_vd = tmp(seq(0,tmp.size()-2));

    cout << "s:  " << m_vs.transpose() << endl;
    cout << "a:  " << m_va.transpose()*180/M_PI << endl;
    cout << "k:  " << m_vk.transpose() << endl;
    cout << "d:  " << m_vd.transpose()*180/M_PI << endl;

    //add 90 deg if not present
    bool is90 = false;
    for(int i=0;i<m_vd.size();i++)
    {
        if (m_vd(i) == Tools::deg2rad(90))
        {
            is90 = true;
            break;
        }
    }
    if (!is90)
    {
        m_vd.conservativeResize(m_vd.size() + 1);
        m_vd(m_vd.size() - 1) = Tools::deg2rad(90);
    }

    MatrixXd uniformSampling(m_vs.size()*
                             5*
                             m_va.size()*
                             m_vk.size()*
                             m_vd.size(),5);

    cout << "Generating: " << m_vs.size()*
                             5*
                             m_va.size()*
                             m_vk.size()*
                             m_vd.size() << " samples" << endl;
    MatrixXd allOrthoResp(m_vs.size()*m_va.size(),4);

    for (int i=0;i<m_vs.size();i++)
    {
        for(int j=0;j<m_va.size();j++)
        {
            allOrthoResp.row(i*m_va.size()+j) << m_vs(i),m_va(j),0,0;
        }
    }
    vector<Vector2d> vorthoResp(m_vs.size()*m_va.size());
    progressbar bar(m_vs.size()*m_va.size());
    for (int i=0;i<vorthoResp.size();i++)
    {
        vorthoResp[i] = orthogonalResponse(allOrthoResp(i,0),allOrthoResp(i,1));
            bar.update();
    }

    int cpt=0;
    for (int i=0;i<m_vs.size();i++)
    {
        for(int j=0;j<m_va.size();j++)
        {
            Vector2d orthResp = vorthoResp[i*m_va.size()+j];
            allOrthoResp.row(i*m_va.size()+j) << m_vs(i),m_va(j),orthResp(0),orthResp(1);

            VectorXd vc = VectorXd::Zero(5);
            vc << orthResp(0)-0.1, 
                  orthResp(0)-0.05,
                  orthResp(0), 
                  orthResp(0)+0.05, 
                  orthResp(0)+0.1;

            for (int k=0;k<vc.size();k++)
            {
                for (int l=0;l<m_vk.size();l++)
                {
                    for (int m=0;m<m_vd.size();m++)
                    {
                        uniformSampling.row(cpt) << m_vs(i),vc(k),m_va(j),m_vk(l),m_vd(m);
                        cpt++;
                    }
                }
            }
        }
    }
    Tools::writeMatrix(allOrthoResp,m_matdir+"/orthResp.csv");
    m_scakd = Tools::deleteDuplicateRows(uniformSampling);
    cout << "Uniform Sampling: " << m_scakd.rows() << endl;
}

Matrix<double,1,6> DataGeneration::derivativesSafe(double s,double c,double a,double k, double d, double h, double &errorConvergence) const
{
    Matrix<double,1,6> FW;

    if( (s==0) and (c==0) and (k==0) )
    {
        double energy = m_slv.getRestEnergy();
        FW << 0,0,0,0,0,energy;
        return FW;
    }
    Solver slv = m_slv;
    double dtemp = (k==0) ? 0 : d;
    MatrixXd Xf = slv.relax(s,c,a,k,dtemp,errorConvergence);
    double energy = slv.getElasticEnergy(Xf,s,c,a,k,d);
    FW(5) = energy;
    double energy_s = slv.getElasticEnergy(Xf,s+h,c,a,k,d);
    FW(0) = (energy_s-energy)/h;
    double energy_c = slv.getElasticEnergy(Xf,s,c+h,a,k,d);
    FW(1) = (energy_c-energy)/h;
    double energy_a = slv.getElasticEnergy(Xf,s,c,a+h,k,d);
    FW(2) = (energy_a-energy)/h;
    double energy_k = slv.getElasticEnergy(Xf,s,c,a,k+h,d);
    FW(3) = (energy_k-energy)/h;
    double energy_d = slv.getElasticEnergy(Xf,s,c,a,k,d+h);
    FW(4) = (energy_d-energy)/h;
    return FW;
}

void DataGeneration::genDerivatives() const
{
    cout << "Generating data" << endl;
    cout << "Writing in:   " << m_matdir+"/mat"+m_matid+"derivatives.csv" << endl;
    if (m_scakd.rows() == 0)
    {
        cout << "DataGeneration::gen::m_scakd is empty" << endl;
        exit(0);
    }
    vector<Matrix<double,1,6>> v_res(m_scakd.rows());
    vector<double> v_error(m_scakd.rows());
    progressbar bar(m_scakd.rows());
    int cpt = 0;
    int cptInvalid = 0;
    #pragma omp parallel for
    for(int i=0;i<m_scakd.rows();i++)
    {
        double error;
        v_res[i] = derivativesSafe(m_scakd(i,0),m_scakd(i,1),m_scakd(i,2),m_scakd(i,3),m_scakd(i,4),1e-6,error);
        v_error[i] = error;
         #pragma omp critical
             bar.update();
    }

    MatrixXd verrors(m_scakd.rows(),1);
    MatrixXd vder(m_scakd.rows(),6);
    MatrixXd vinput(m_scakd.rows(),5);

    for(int i=0;i<m_scakd.rows();i++)
    {
        vder.row(i) = v_res[i];
        vinput.row(i) = m_scakd.row(i);
        verrors(i) = v_error[i];
    }

    //divide by volume
    Solver slv = m_slv;
    double volume = slv.getVolume();
    cout << "Volume : " << volume << endl;
    if (isnan(vder.norm()))
    {
        cout << "vder is nan" << endl;
        exit(0);
    }

    //write
    Tools::writeMatrix(vder/volume,m_matdir+"/mat"+m_matid+"derivatives.csv");
    Tools::writeMatrix(vinput,m_matdir+"/mat"+m_matid+"scakd.csv");
    Tools::writeMatrix(verrors,m_matdir+"/mat"+m_matid+"errors.csv");
    cout << "MaxError: " << verrors.maxCoeff() << endl;
}

Vector2d DataGeneration::findMinEnergy(double s,double &c, double a,double &cmin,double &cmax) const
{

    unsigned int nthreads = 12;
    ArrayXd vc = ArrayXd::LinSpaced(nthreads, cmin, cmax);
    double h = abs(vc(0)-vc(1));
    vector<double> energies(nthreads);
    vector<Solver> v_slv(nthreads);
    for(int i=0;i<nthreads;i++)
        {v_slv[i] = m_slv;}
    #pragma omp parallel for
    for (int i = 0; i < vc.size(); i++)
    {
        bool success;
        double energy = v_slv[i].getMinimalEnergy(s,vc(i),a,0,0,success);
        if (success)
            {energies[i] = energy;}
        else
            {energies[i] = 1e12;}
    }

    int minIdx = std::min_element(energies.begin(), energies.end()) - energies.begin();
    if (energies[minIdx] == 1e12)
    {
        cout << "parallelSweep::fail" << endl;
        exit(0);
    }
    cmin = vc(minIdx)-h;
    cmax = vc(minIdx)+h;
    c = vc(minIdx);
    double energy = energies[minIdx];
    double gradient;
    if (minIdx == 0)
    {
        gradient = abs(energies[minIdx + 1] - energies[minIdx]) / abs(h);
    }
    else if (minIdx==vc.size()-1)
    {
        gradient = abs(energies[minIdx - 1] - energies[minIdx]) / abs(h);
    }
    else
    {
        gradient = abs(min(energies[minIdx - 1], energies[minIdx + 1]) - energies[minIdx]) / abs(h);
    }
    Vector2d res;
    res << energy, gradient;
    return res;
}

Vector2d DataGeneration::orthogonalResponse(double s,double a) const
{
    int nbMaxReps = 25;
    VectorXd vorthSearch(nbMaxReps);
    VectorXd vgradient(nbMaxReps);
    VectorXd vorthoEnergies(nbMaxReps);
    vorthSearch.setOnes(); vorthSearch *= 1e12;
    vorthoEnergies.setOnes(); vorthoEnergies *= 1e12;
    double cmin = -0.8;
    double cmax = 0.8;
    double c = 1e10;
    for(int k=0;k<nbMaxReps;k++)
    {
        Vector2d EG = findMinEnergy(s,c, a,cmin,cmax);
        vorthSearch(k) = c;
        vorthoEnergies(k) = EG(0);
        vgradient(k) = EG(1);
        if (EG(1) < 1.0)
        {
            break;
        }
    }
    Index minIdx;
    vorthoEnergies.minCoeff(&minIdx);
    Vector2d res;
    res << vorthSearch(minIdx), vgradient(minIdx);
    return res;
}

void DataGeneration::orthRespProfile(double s,double a) const
{
    int nbSamples = 100;
    VectorXd vc = ArrayXd::LinSpaced(nbSamples, -0.8, 0.8);
    VectorXd vorthResp(vc.size());
    vector<Solver> v_slv(nbSamples);
    for(int i=0;i<nbSamples;i++)
        {v_slv[i] = m_slv;}
    #pragma omp parallel for
    for(int i=0;i<vc.size();i++)
    {
        bool success;
        vorthResp(i) =  v_slv[i].getMinimalEnergy(s,vc(i),a,0,0,success);
        if (!success)
        {
            vorthResp(i) = 1e12;
        }
    }
    MatrixXd res(vc.size(),2);
    res.col(0) = vc;
    res.col(1) = vorthResp;
    Tools::writeMatrix(res,"orthResp.csv");
}
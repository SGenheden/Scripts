#include <iostream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

void calc_zero(double r, double &fx, double &fdx ) {
    fx = 0.0;
    fdx = 0.0;
}

void calc_zerolj(double r, double &gx, double &gdx, double &hx, double &hdx  ) {
    gx = 0.0;
    gdx = 0.0;
    hx = 0.0;
    hdx = 0.0;
}

void calc_coul(double r, double dielectric, double &fx, double &fdx ) {
    fx = 1.0 / r;
    fdx = 1.0 / (r*r);
    fx /= dielectric;
    fdx /= dielectric;
}

void calc_rf(double r, double cutoff, double dielectric, double epsilon, double &fx, double &fdx ) {

    double one = 1.0;
    double two = 2.0;
    double three = 3.0;

    double cutoff2 = cutoff*cutoff;
    double cutoff3 = cutoff2*cutoff;
    double krf = (epsilon-dielectric)/(cutoff3*(two*epsilon+dielectric));
    double crf = (three*epsilon)/(cutoff*(two*epsilon+dielectric));
    double fshift = (three*dielectric)/(cutoff*cutoff*(two*epsilon+dielectric));

    double rinv = one / r;
    double rinv2 = rinv*rinv;
    double r2 = r * r;

    fx = rinv + krf*r2 - crf;
    fdx = rinv2 - two*krf*r - fshift;

    fx /= dielectric;
    fdx /= dielectric;
}

void calc_cutlj(double r, double &gx, double &gdx, double &hx, double &hdx) {
    double rinv = 1.0 / r;
    double rinv2 = rinv*rinv;
    double rinv6 = rinv2*rinv2*rinv2;
    double rinv12 = rinv6*rinv6;
    double rinv7  = rinv6*rinv;
    double rinv13 = rinv12*rinv;

    gx = -rinv6;
    gdx = -6.0*rinv7;

    hx = rinv12;
    hdx = 12.0*rinv13;
}

void calc_shiftedlj(double r, double cutoff, double &gx, double &gdx, double &hx, double &hdx) {

    double zero = 0.0;
    double one = 1.0;
    double two = 2.0;
    double three = 3.0;
    double third = one/three;
    double four = 4.0;
    double fourth = one/four;
    double six = 6.0;
    double twelve = 12.0;

    double power1 = one;
    double power2 = six;
    double power3 = twelve;

    double shiftlj = 0.90;
    double rmodlj = cutoff - shiftlj;
    double rmodlj2 = rmodlj*rmodlj;
    double rmodlj3 = rmodlj2*rmodlj;
    double rmodlj4 = rmodlj2*rmodlj2;

    double a2 = -((power2+four)*cutoff - (power2+one)*shiftlj)/
        ((pow(cutoff,power2+two))*rmodlj2) ;
    double b2 = ((power2+three)*cutoff - (power2+one)*shiftlj)/
        ((pow(cutoff,power2+two))*rmodlj3);
    double c2 = one/pow(cutoff,power2) - third*power2*a2*rmodlj3
        - fourth*power2*b2*rmodlj4;
    double a3 = -((power3+four)*cutoff - (power3+one)*shiftlj)/
        ((pow(cutoff,power3+two))*rmodlj2);
    double b3 = ((power3+three)*cutoff - (power3+one)*shiftlj)/
        ((pow(cutoff,power3+two))*rmodlj3);
    double c3 = one/pow(cutoff,power3) - third*power3*a3*rmodlj3
        - fourth*power3*b3*rmodlj4;

    double rlj = r - shiftlj;
    double rlj2 = rlj*rlj;
    double rlj3 = rlj2*rlj;
    double rlj4 = rlj2*rlj2;

    if (r <= shiftlj) {
        gx = -one/pow(r,power2) + c2;
        gdx = -power2/pow(r,power2+one);
        hx = one/pow(r,power3) - c3;
        hdx = power3/pow(r,power3+one);
    } else {
        gx = -one/pow(r,power2) + third*power2*a2*rlj3
            + fourth*power2*b2*rlj4 + c2;
        gdx = -power2/pow(r,power2+one) - power2*(a2*rlj2 + b2*rlj3);
        hx = one/pow(r,power3) - third*power3*a3*rlj3
            - fourth*power3*b3*rlj4 - c3;
        hdx = power3/pow(r,power3+one) + power3*(a3*rlj2 + b3*rlj3);
    }
}

int main()
{

    std::cout << "Enter the type of potential for electrostatics: coulomb (coul) or reaction-field (rf):";
    string ele = "coul";
    cin >> ele;

    std::cout << "Enter the type of potential for van der Waals: cut-off (cut) or shifted (shifted): ";
    string vdw = "cut";
    cin >> vdw;

    std::cout << "Enter the cutoff: ";
    double cutoff = 1.200;
    cin >> cutoff;

    std::cout << "Enter the dielectric constant: ";
    double dielectric = 1.0;
    cin >> dielectric;

    double epsilon = 1.0;
    if (ele == "rf") {
        std::cout << "Enter the RF epsilon: ";
        cin >> epsilon;
    }

    std::cout << "Enter the output filename: ";
    string outname = "";
    cin >> outname;

    cout << "\n";

    double delta = 0.002000;
    FILE *fout = fopen(outname.c_str(), "w");
    for (double r = 0.0; r<cutoff+2.0*delta+1.0; r+=delta) {
        if ((r < 4.0E-2) || (r > cutoff+delta)) {
            fprintf(fout, "%12.10e   %12.10e %12.10e   %12.10e %12.10e   %12.10e %12.10e\n", r,0.0,0.0,0.0,0.0,0.0,0.0);
        } else {
            double fx, fdx, gx, gdx, hx, hdx = 0.0;
            // Calculate electrostatics
            if (ele == "coul") {
                calc_coul(r, dielectric, fx, fdx);
            } else if (ele == "rf") {
                calc_rf(r, cutoff, dielectric, epsilon, fx, fdx);
            } else {
                calc_zero(r, fx, fdx);
            }
            // Calculate van der Waals potential
            if (vdw == "cut") {
                calc_cutlj(r, gx, gdx, hx, hdx);
            } else if (vdw == "shifted") {
                calc_shiftedlj(r, cutoff, gx, gdx, hx, hdx);
            } else   {
                calc_zerolj(r, gx, gdx, hx, hdx);
            }
            fprintf(fout, "%12.10e   %12.10e %12.10e   %12.10e %12.10e   %12.10e %12.10e\n", r,fx, fdx, gx, gdx, hx, hdx);
        }
    }

    return 0;

}

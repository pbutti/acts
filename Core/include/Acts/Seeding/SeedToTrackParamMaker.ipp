#include "Acts/Seeding/SeedToTrackParamMaker.hpp"

#include <cmath>


namespace Acts {
    
    /// V. Karimaki NIM A305 (1991) 187-191 - no weights are used
    /// In Karimaki's fit, d0 is the distance of the closest approach to the origin, 1/R is the curvature, 
    /// phi is the angle of the direction propagation (counter clock wise as positive)
    
    template <typename external_spacepoint_t>
    bool SeedToTrackParamMaker::KarimakiFit(const std::vector<external_spacepoint_t*>&sp, std::array<double,9>& data) {
        
        if (!sp.size())
            return false;
        
        double x2m=0.,  xm=0.   ;
        double xym=0.           ;
        double y2m=0.,  ym=0.   ;
        double r2m=0.,  r4m=0.  ;
        double xr2m=0., yr2m=0. ;

        int n = sp.size();

        for (size_t i_sp = 0; i_sp<sp.size(); i_sp++) {
            
            double x  = sp[i_sp]->x();
            double y  = sp[i_sp]->y();
            double r2 = x*x+y*y;
            
            x2m  += x*x   * 1./n;
            xm   += x     * 1./n;
            xym  += x*y   * 1./n;
            y2m  += y*y   * 1./n;
            ym   += y     * 1./n;
            r2m  += r2    * 1./n;
            r4m  += r2*r2 * 1./n;
            xr2m += x*r2  * 1./n;
            yr2m += y*r2  * 1./n;
        }
        
        double Cxx   = x2m  - xm*xm;
        double Cxy   = xym  - xm*ym;
        double Cyy   = y2m  - ym*ym;
        double Cxr2  = xr2m - xm*r2m; 
        double Cyr2  = yr2m - ym*r2m;
        double Cr2r2 = r4m  - r2m*r2m;
        
        double q1 = Cr2r2*Cxy - Cxr2*Cyr2 ;
        double q2 = Cr2r2*(Cxx-Cyy) - Cxr2*Cxr2+Cyr2*Cyr2;
        
        double phi   = 0.5 * atan(2*q1/q2);
        double k     = (sin(phi)*Cxr2 - cos(phi)*Cyr2) * (1./Cr2r2);
        double delta = -k*r2m + sin(phi)*xm - cos(phi)*ym;
        
        double rho = (2*k) / (sqrt(1-4*delta*k));
        double d   = (2*delta) / (1+ sqrt(1-4*delta*k));
        
        std::cout<<"rho="<<rho<<" d="<<d<<" phi="<<phi<<std::endl; 
        
        return true;
        
    }
    
    /// see https://acode-browser.usatlas.bnl.gov/lxr/source/athena/InnerDetector/InDetRecTools/SiTrackMakerTool_xk/src/SiTrackMaker_xk.cxx
    template <typename external_spacepoint_t>
    bool SeedToTrackParamMaker::FitSeedAtlas(const Seed<external_spacepoint_t>& seed, std::array<double,9>& data, const Transform3D& Tp, const double& bFieldZ) {
        return FitSeedAtlas(seed.sp(),data, Tp, bFieldZ);
    }

        
    //double H = 0.0015; //kTesla 

    /// This method gives an estimate of the track parameters of a seed using a conformal map transformation
    /// The track parameters are of the form l1, l2, phi, theta, q/p or d0, z0, phi, theta, q/p. 
    /// phi0 is the angle of the track direction with respect the origin, positive when counter clock-wise
    
    
    template <typename external_spacepoint_t>
    bool SeedToTrackParamMaker::FitSeedAtlas(const std::vector<external_spacepoint_t>& sp, std::array<double,9>& data, const Transform3D& Tp, const double& bFieldZ) {
        
        m_pTmin = 0.1;
        
        /// Define the locations of the space points with respect to the first one in the Space point vector
        
        double x0 = sp[0]->x();
        double y0 = sp[0]->y();
        double z0 = sp[0]->z();

        double r0 = sqrt(x0*x0+y0*y0);
        
        double x1 = sp[1]->x() - x0;
        double y1 = sp[1]->y() - y0;
        double x2 = sp[2]->x() - x0;
        double y2 = sp[2]->y() - y0;
        double z2 = sp[2]->z() - z0; 
        

        /// Define conformal map variables
        double u1 = 1./sqrt(x1*x1+y1*y1)       ;
        double rn = x2*x2+y2*y2                ;
        double r2 = 1./rn                      ;
        
        /// a = cos(phi_0), b = sin(phi_0) in some other notations
        /// and solve
        double a  = x1*u1                      ;
        double b  = y1*u1                      ;
        double u2 = (a*x2+b*y2)*r2             ;
        double v2 = (a*y2-b*x2)*r2             ;
        double A  = v2/(u2-u1)                 ;
        double B  = 2.*(v2-A*u2)               ;
        
        /// 1/helixradius = C
        double C  = B/sqrt(1.+A*A)             ;
                
        double T  = z2*sqrt(r2)/(1.+.04*C*C*rn);
        
        /// Project to the surface
        double Ax[3] = {Tp(0,0),Tp(1,0),Tp(2,0)};
        double Ay[3] = {Tp(0,1),Tp(1,1),Tp(2,1)};
        double D [3] = {Tp(0,3),Tp(1,3),Tp(2,3)};
        
        double d[3]  = {x0-D[0], y0-D[1], z0-D[2]};
        
        /// l1 is the (most) sensitive direction, l2 is the un- (less) sensitive direction
        data[0] = d[0]*Ax[0]+d[1]*Ax[1]+d[2]*Ax[2]; 
        data[1] = d[0]*Ay[0]+d[1]*Ay[1]+d[2]*Ay[2]; 
        
        data[2] = std::atan2(b+a*A,a-b*A);
        data[3] = std::atan2(1.,T);
        data[5] = -C / (300. * bFieldZ);
        
        /// B in this computation is twice the usual B
        double b_c = 1/B;
        double a_c = -1*b_c*A;
        double R = 1./C;
        
        std::cout<<"(a_c, b_c) = "<<a_c<<","<<b_c<<std::endl;
        std::cout<<"R="<<R<<std::endl;
        
        // wrt real origin
        double a_cp = a_c*(a) - b_c*(b) + x0;
        double b_cp = b_c*(a) + a_c*(b) + y0;
        double ip = (1./(2.*R))*(a_cp*a_cp+b_cp*b_cp -R*R);
        
        std::cout<<"SeedToTrackParamMaker (a_cp,b_cp)="<<a_cp<<","<<b_cp<<std::endl;
        std::cout<<"Ip="<<ip<<std::endl;
        
        /// Return false if the transverse momentum is less than 90% the minimum momentum
        if (fabs(data[5]*m_pTmin) > 1.1)
            return false;
        
        /// Momentum is given by pT/cos(theta)
        data[4] = data[5] / std::sqrt(1. + T*T);

        /// Store the reference point
        data[6] = x0;
        data[7] = y0;
        data[8] = z0;
        
        return true;
    }
    
    template <typename external_spacepoint_t>
    bool SeedToTrackParamMaker::FitSeedLinCircle(const Seed<external_spacepoint_t>& seed, std::vector<double>& data) {
        return true;
    }
    
    template <typename external_spacepoint_t>
    bool SeedToTrackParamMaker::FitSeedLinPar (const Seed<external_spacepoint_t>& seed, std::vector<double>& data)  {
        return true;
    }
    
    template <typename external_spacepoint_t>
    bool SeedToTrackParamMaker::LinearParabolaFit(std::vector<external_spacepoint_t>& sPs) {
        return true;
    }
        


}

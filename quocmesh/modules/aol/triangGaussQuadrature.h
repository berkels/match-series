#ifndef __TRIANGGAUSSQUADRATURE_H
#define __TRIANGGAUSSQUADRATURE_H

#include <smallVec.h>

namespace aol {


//! Different quadrature types for triangular meshes


/*!
 * \author Droske
 */
// CenterQuadrature
template <typename RealType>
class CenterQuadrature {
public:
  enum { numQuadPoints = 1 };

  inline const aol::Vec2<RealType>& getRefCoord ( int ) const  {
    static const aol::Vec2<RealType> c ( 1. / 3., 1. / 3. );
    return c;
  }

  inline RealType getWeight ( int ) const {
    return 1.;
  }
};

/*!
 * \author Droske
 */
// TriQuadrature
template <typename RealType>
class TriQuadrature {
public:
  enum { numQuadPoints = 3 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[3] = {
                                              aol::Vec2<RealType> ( 1. / 6., 1. / 6. ),
                                              aol::Vec2<RealType> ( 1. / 6., 2. / 3. ),
                                              aol::Vec2<RealType> ( 2. / 3., 1. / 6. )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int ) const {
    return 1. / 3.;
  }
};


/*!
 * \author Perl
 */
// EdgeQuadrature
template <typename RealType>
class EdgeQuadrature {
public:
  enum { numQuadPoints = 3 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[3] = {
                                              aol::Vec2<RealType> ( 1./2., 0. ),
                                              aol::Vec2<RealType> ( 0., 1./2. ),
                                              aol::Vec2<RealType> ( 1./2., 1./2. )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int ) const {
    return 1./3.;
  }
};

/*!
 * \author Perl
 */
// EdgeQuadrature
template <typename RealType>
class CornerQuadrature {
public:
  enum { numQuadPoints = 3 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[3] = {
                                              aol::Vec2<RealType> ( 0., 0. ),
                                              aol::Vec2<RealType> ( 0., 1. ),
                                              aol::Vec2<RealType> ( 1., 0. )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int ) const {
    return 1./3.;
  }
};

/*!
 * \author Perl
 */
// QuadriQuadrature, exact for polynomials of degree 3
template <typename RealType>
class QuadriQuadrature {
public:
  enum { numQuadPoints = 4 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[4] = {
                                              aol::Vec2<RealType> ( 1. / 3., 1. / 3. ),
                                              aol::Vec2<RealType> ( 3. / 5., 1. / 5. ),
                                              aol::Vec2<RealType> ( 1. / 5., 1. / 5. ),
                                              aol::Vec2<RealType> ( 1. / 5., 3. / 5. )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[4] = { 
                    RealType (-27./48. ),
                    RealType ( 25./48. ),
                    RealType ( 25./48. ),
                    RealType ( 25./48. )
                  };
    return w[quadpoint];
  }
};

/*!
 * \author Perl
 */
// SexaQuadrature, exact for polynomials of degree 4
template <typename RealType>
class SexaQuadrature {
public:
  enum { numQuadPoints = 6 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[6] = {
                                              aol::Vec2<RealType> ( 0.44594849091597 , 0.44594849091597 ),
                                              aol::Vec2<RealType> ( 0.44594849091597 , 0.10810301816807 ),
                                              aol::Vec2<RealType> ( 0.10810301816807 , 0.44594849091597 ),
                                              aol::Vec2<RealType> ( 0.09157621350977 , 0.09157621350977 ),
                                              aol::Vec2<RealType> ( 0.09157621350977 , 0.81684757298046 ),
                                              aol::Vec2<RealType> ( 0.81684757298046 , 0.09157621350977 )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[6] = { 
                    RealType ( 0.22338158967801 ),
                    RealType ( 0.22338158967801 ),
                    RealType ( 0.22338158967801 ),
                    RealType ( 0.10995174365532 ),
                    RealType ( 0.10995174365532 ),
                    RealType ( 0.10995174365532 )
                  };
    return w[quadpoint];
  }
};

/*!
 * \author Perl
 */
// SeptiQuadrature, exact for polynomials of degree 5
template <typename RealType>
class SeptiQuadrature {
public:
  enum { numQuadPoints = 7 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[7] = {
                                              aol::Vec2<RealType> ( 0.33333333333333 , 0.33333333333333 ),
                                              aol::Vec2<RealType> ( 0.47014206410511 , 0.47014206410511 ),
                                              aol::Vec2<RealType> ( 0.47014206410511 , 0.05971587178977 ),
                                              aol::Vec2<RealType> ( 0.05971587178977 , 0.4701420641051 ),
                                              aol::Vec2<RealType> ( 0.10128650732346 , 0.10128650732346 ),
                                              aol::Vec2<RealType> ( 0.10128650732346 , 0.79742698535309 ),
                                              aol::Vec2<RealType> ( 0.79742698535309 , 0.10128650732346 )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[7] = { 
                    RealType ( 0.22500000000000 ),
                    RealType ( 0.13239415278851 ),
                    RealType ( 0.13239415278851 ),
                    RealType ( 0.13239415278851 ),
                    RealType ( 0.12593918054483 ),
                    RealType ( 0.12593918054483 ),
                    RealType ( 0.12593918054483 )
                  };
    return w[quadpoint];
  }
};

/*!
 * \author Perl
 */
// DuodecQuadrature, exact for polynomials of degree 6
template <typename RealType>
class DuodecQuadrature {
public:
  enum { numQuadPoints = 12 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[12] = {
                                              aol::Vec2<RealType> ( 0.24928674517091 , 0.24928674517091 ),
                                              aol::Vec2<RealType> ( 0.24928674517091 , 0.50142650965818 ),
                                              aol::Vec2<RealType> ( 0.50142650965818 , 0.24928674517091 ),
                                              aol::Vec2<RealType> ( 0.06308901449150 , 0.06308901449150 ),
                                              aol::Vec2<RealType> ( 0.06308901449150 , 0.87382197101700 ),
                                              aol::Vec2<RealType> ( 0.87382197101700 , 0.06308901449150 ),
                                              aol::Vec2<RealType> ( 0.31035245103378 , 0.63650249912140 ),
                                              aol::Vec2<RealType> ( 0.63650249912140 , 0.05314504984482 ),
                                              aol::Vec2<RealType> ( 0.05314504984482 , 0.31035245103378 ),
                                              aol::Vec2<RealType> ( 0.63650249912140 , 0.31035245103378 ),
                                              aol::Vec2<RealType> ( 0.31035245103378 , 0.05314504984482 ),
                                              aol::Vec2<RealType> ( 0.05314504984482 , 0.63650249912140 )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[12] = { 
                    RealType ( 0.11678627572638 ),
                    RealType ( 0.11678627572638 ),
                    RealType ( 0.11678627572638 ),
                    RealType ( 0.05084490637021 ),
                    RealType ( 0.05084490637021 ),
                    RealType ( 0.05084490637021 ),
                    RealType ( 0.08285107561837 ),
                    RealType ( 0.08285107561837 ),
                    RealType ( 0.08285107561837 ),
                    RealType ( 0.08285107561837 ),
                    RealType ( 0.08285107561837 ),
                    RealType ( 0.08285107561837 )
                  };
    return w[quadpoint];
  }
};

/*!
 * \author Perl
 */
// TredecQuadrature, exact for polynomials of degree 6
template <typename RealType>
class TredecQuadrature {
public:
  enum { numQuadPoints = 13 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[13] = {
                                              aol::Vec2<RealType> ( 0.33333333333333 , 0.33333333333333 ),
                                              aol::Vec2<RealType> ( 0.26034596607904 , 0.26034596607904 ),
                                              aol::Vec2<RealType> ( 0.26034596607904 , 0.47930806784192 ),
                                              aol::Vec2<RealType> ( 0.47930806784192 , 0.26034596607904 ),
                                              aol::Vec2<RealType> ( 0.06513010290222 , 0.06513010290222 ),
                                              aol::Vec2<RealType> ( 0.06513010290222 , 0.86973979419557 ),
                                              aol::Vec2<RealType> ( 0.86973979419557 , 0.06513010290222 ),
                                              aol::Vec2<RealType> ( 0.31286549600487 , 0.63844418856981 ),
                                              aol::Vec2<RealType> ( 0.63844418856981 , 0.04869031542532 ),
                                              aol::Vec2<RealType> ( 0.04869031542532 , 0.31286549600487 ),
                                              aol::Vec2<RealType> ( 0.63844418856981 , 0.31286549600487 ),
                                              aol::Vec2<RealType> ( 0.31286549600487 , 0.04869031542532 ),
                                              aol::Vec2<RealType> ( 0.04869031542532 , 0.63844418856981 )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[13] = { 
                    RealType ( -0.14957004446768 ),
                    RealType ( 0.17561525743321 ),
                    RealType ( 0.17561525743321 ),
                    RealType ( 0.17561525743321 ),
                    RealType ( 0.05334723560884 ),
                    RealType ( 0.05334723560884 ),
                    RealType ( 0.05334723560884 ),
                    RealType ( 0.07711376089026 ),
                    RealType ( 0.07711376089026 ),
                    RealType ( 0.07711376089026 ),
                    RealType ( 0.07711376089026 ),
                    RealType ( 0.07711376089026 ),
                    RealType ( 0.07711376089026 )
                  };
    return w[quadpoint];
  }
};

/*!
 * \author Perl
 */
// SexdecQuadrature, exact for polynomials of degree 8
template <typename RealType>
class SexdecQuadrature {
public:
  enum { numQuadPoints = 16 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[16] = {
    aol::Vec2<RealType> ( 0.333333333333333333333333333333333 , 0.333333333333333333333333333333333 ),
    aol::Vec2<RealType> ( 0.459292588292723156028815514494169 , 0.459292588292723156028815514494169 ),
    aol::Vec2<RealType> ( 0.459292588292723156028815514494169 , 0.081414823414553687942368971011661 ),
    aol::Vec2<RealType> ( 0.081414823414553687942368971011661 , 0.459292588292723156028815514494169 ),
    aol::Vec2<RealType> ( 0.170569307751760206622293501491464 , 0.170569307751760206622293501491464 ),
    aol::Vec2<RealType> ( 0.170569307751760206622293501491464 , 0.658861384496479586755412997017071 ),
    aol::Vec2<RealType> ( 0.658861384496479586755412997017071 , 0.170569307751760206622293501491464 ),
    aol::Vec2<RealType> ( 0.050547228317030975458423550596598 , 0.050547228317030975458423550596598 ),
    aol::Vec2<RealType> ( 0.050547228317030975458423550596598 , 0.898905543365938049083152898806802 ),
    aol::Vec2<RealType> ( 0.898905543365938049083152898806802 , 0.050547228317030975458423550596598 ),
    aol::Vec2<RealType> ( 0.263112829634638113421785786284643 , 0.728492392955404281241000379176062 ),
    aol::Vec2<RealType> ( 0.728492392955404281241000379176062 , 0.008394777409957605337213834539294 ),
    aol::Vec2<RealType> ( 0.008394777409957605337213834539294 , 0.263112829634638113421785786284643 ),
    aol::Vec2<RealType> ( 0.728492392955404281241000379176062 , 0.263112829634638113421785786284643 ),
    aol::Vec2<RealType> ( 0.263112829634638113421785786284643 , 0.008394777409957605337213834539294 ),
    aol::Vec2<RealType> ( 0.008394777409957605337213834539294 , 0.728492392955404281241000379176062 )
    };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[16] = { 
                    RealType ( 0.144315607677787168251091110489064 ),
                    RealType ( 0.095091634267284624793896104388584 ),
                    RealType ( 0.095091634267284624793896104388584 ),
                    RealType ( 0.095091634267284624793896104388584 ),
                    RealType ( 0.103217370534718250281791550292129 ),
                    RealType ( 0.103217370534718250281791550292129 ),
                    RealType ( 0.103217370534718250281791550292129 ),
                    RealType ( 0.032458497623198080310925928341780 ),
                    RealType ( 0.032458497623198080310925928341780 ),
                    RealType ( 0.032458497623198080310925928341780 ),
                    RealType ( 0.027230314174434994264844690073908 ),
                    RealType ( 0.027230314174434994264844690073908 ),
                    RealType ( 0.027230314174434994264844690073908 ),
                    RealType ( 0.027230314174434994264844690073908 ),
                    RealType ( 0.027230314174434994264844690073908 ),
                    RealType ( 0.027230314174434994264844690073908 )
                  };
    return w[quadpoint];
  }
};

/*!
 * \author Perl
 */
// GaussDegree8Quadrature, exact for polynomials of degree 10
template <typename RealType>
class GaussDegree10Quadrature {
public:
  enum { numQuadPoints = 25 };

  inline const aol::Vec2<RealType>& getRefCoord ( int quadpoint ) const  {
    static const aol::Vec2<RealType> c[25] = {
    aol::Vec2<RealType> ( 0.333333333333333 , 0.333333333333333 ),
    aol::Vec2<RealType> ( 0.028844733232685 , 0.485577633383657 ),
    aol::Vec2<RealType> ( 0.485577633383657 , 0.028844733232685 ),
    aol::Vec2<RealType> ( 0.485577633383657 , 0.485577633383657 ),
    aol::Vec2<RealType> ( 0.781036849029926 , 0.109481575485037 ),
    aol::Vec2<RealType> ( 0.109481575485037 , 0.781036849029926 ),
    aol::Vec2<RealType> ( 0.109481575485037 , 0.109481575485037 ),
    aol::Vec2<RealType> ( 0.141707219414880 , 0.307939838764121 ),
    aol::Vec2<RealType> ( 0.141707219414880 , 0.550352941820999 ),
    aol::Vec2<RealType> ( 0.307939838764121 , 0.141707219414880 ),
    aol::Vec2<RealType> ( 0.307939838764121 , 0.550352941820999 ),
    aol::Vec2<RealType> ( 0.550352941820999 , 0.141707219414880 ),
    aol::Vec2<RealType> ( 0.550352941820999 , 0.307939838764121 ),
    aol::Vec2<RealType> ( 0.025003534762686 , 0.246672560639903 ),
    aol::Vec2<RealType> ( 0.025003534762686 , 0.728323904597411 ),
    aol::Vec2<RealType> ( 0.246672560639903 , 0.025003534762686 ),
    aol::Vec2<RealType> ( 0.246672560639903 , 0.728323904597411 ),
    aol::Vec2<RealType> ( 0.728323904597411 , 0.025003534762686 ),
    aol::Vec2<RealType> ( 0.728323904597411 , 0.246672560639903 ),
    aol::Vec2<RealType> ( 0.009540815400299 , 0.066803251012200 ),
    aol::Vec2<RealType> ( 0.009540815400299 , 0.923655933587500 ),
    aol::Vec2<RealType> ( 0.066803251012200 , 0.009540815400299 ),
    aol::Vec2<RealType> ( 0.066803251012200 , 0.923655933587500 ),
    aol::Vec2<RealType> ( 0.923655933587500 , 0.009540815400299 ),
    aol::Vec2<RealType> ( 0.923655933587500 , 0.066803251012200 )
    };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[25] = { 
                    RealType ( 0.090817990382754 ),
                    RealType ( 0.036725957756467 ),
                    RealType ( 0.036725957756467 ),
                    RealType ( 0.036725957756467 ),
                    RealType ( 0.045321059435528 ),
                    RealType ( 0.045321059435528 ),
                    RealType ( 0.045321059435528 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 )
                  };
    return w[quadpoint];
  }
};

}

#endif

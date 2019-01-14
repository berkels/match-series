#ifndef __INTERRUPTABLEITERATIVEINFO_H
#define __INTERRUPTABLEITERATIVEINFO_H

#include <ctrlCCatcher.h>
#include <iterativeInfo.h>
#include <solverInfo.h>
#include <NewtonInfo.h>

namespace aol {

/****************************************************************************
 *
 *        CLASS InterruptableIterativeInfo
 */
/**
 *  \brief Decorator for any IterativeInfo subclass that enables
 *         catching Ctrl-C pressures.
 *
 *  The class is derived from its template argument and thus serves as
 *  decorator (GOF design pattern) for an arbitrary IterativeInfo subclass.
 *  Simply create (e.g.) an InterruptableInterativeInfo<NewtonInfo> and pass
 *  it to your Newton iteration object. Then, when pressing Ctrl-C, your
 *  program will not be terminated, but normally continued. At the end of the
 *  current step, an alert message shows up and asks if you want to terminate
 *  the current iteration. If so, the iteration process will be finished as
 *  if the maximum number of iterations would be exceeded. If you decide not
 *  to terminate your iterations, they will be continued as if nothing would
 *  have happened. While the program waits for your yes/no input to continue
 *  or terminate, "Ctrl-C" has its normal meaning, i. e. you can cancel the
 *  execution immediately.
 *
 *  The concept is most useful when dealing with nested iterations, e.g.
 *  a loop over different grid levels, a Newton iteration for each and
 *  an iteration solver inside Newton's method. When you cancel the solver's
 *  iteration, the Newton scheme will proceed until the end of the current
 *  Newton step. As the iterative solver has done before, it recognized at
 *  the end of its step that Ctrl-C was pressed and asks if you would also
 *  like to end the Newton iterations. This way, your program will correctly
 *  execute for the next grid level, even though the Newton scheme on the
 *  previous level would have taked endless time normally (due to slow
 *  solver's convergence or a large maximum step number).
 *
 *  \author von Deylen, Jan 2009.
 */
template <typename BaseInfo>
class InterruptableIterativeInfo : public BaseInfo {
public:
  typedef typename BaseInfo::RealType RealType;

  //! \short constructor from base info class
  InterruptableIterativeInfo ( const BaseInfo & base );

  //! constructor for BaseInfo = SolverInfo.
  //! will cause compiler errors when called for another BaseInfo type.
  InterruptableIterativeInfo ( RealType TOL,
                               int maxIter,
                               StoppingMode mode = STOPPING_UNSET,
                               bool quietMode = false,
                               ostream & _out = clog,
                               int numAveragingSteps = 10 );

  //! constructor for BaseInfo = NewtonInfo
  //! will cause compiler errors when called for another BaseInfo type.
  InterruptableIterativeInfo ( RealType TOL_Newton,
                               int maxIterNewton,
                               SolverInfo<RealType> & solverInfo,
                               double tauMin = 1E-10,
                               double tauMax = 2.,
                               int numAveragingSteps = 1 );

  //! constructor for BaseInfo = NewtonInfo
  //! will cause compiler errors when called for another BaseInfo type.
  InterruptableIterativeInfo ( RealType TOL_Newton,
                               int maxIterNewton,
                               RealType TOL_solver,
                               int maxIterSolver,
                               StoppingMode mode = STOPPING_UNSET,
                               ostream & out = clog,
                               double tauMin = 1E-10,
                               double tauMax = 2.,
                               int numAveragingSteps = 1 );

  //! installs our own Ctrl-C handler and calls BaseInfos'
  //! startIterations().
  void startIterations();

  //! calls BaseInfos' finishIterations() and removes
  //! our own Ctrl-C handler.
  void finishIterations();

  //! returns true if BaseInfo::maxIterIsExceeded() or
  //! wantsIterrupt() returns true.
  bool maxIterIsExceeded() const;
  //! returns true if BaseInfo::maxIterIsReached() or
  //! wantsIterrupt() returns true.
  bool maxIterIsReached() const;

  //! returns true if "Ctrl-C pressed" flag is risen.
  bool wasInterrupted() const;

protected:
  //! if Ctrl-C was pressed, this functions prints an alert messages
  //! and asks the user to decide if he really wants to terminate
  //! the iteration loop.
  bool wantsInterrupt() const;

  mutable bool _wasInterrupted;
  //! stored the Ctrl-C handler that was used before iterations' start.
  sigfunc _previousCtrlCHandler;
};

/**** IMPLEMENTATION ********************************************************
 *
 *
 *
 */

//---------------------------------------------------------------------------

template <typename BaseInfo>
InterruptableIterativeInfo<BaseInfo>::InterruptableIterativeInfo
              ( const BaseInfo & base )
    : BaseInfo ( base )
    , _wasInterrupted (false)
    , _previousCtrlCHandler ( ErrorSignal )
{}
//---------------------------------------------------------------------------

template <typename BaseInfo>
InterruptableIterativeInfo<BaseInfo>::InterruptableIterativeInfo
              ( RealType TOL, int maxIter,
                StoppingMode mode,
                bool quietMode,
                ostream & out,
                int numAveragingSteps )
    : BaseInfo ( TOL, maxIter, mode, quietMode, out, numAveragingSteps )
    , _wasInterrupted (false)
    , _previousCtrlCHandler ( ErrorSignal )
{}
//---------------------------------------------------------------------------

template <typename BaseInfo>
InterruptableIterativeInfo<BaseInfo>::InterruptableIterativeInfo
              ( RealType TOL_Newton,
                int maxIterNewton,
                SolverInfo<RealType> & solverInfo,
                double tauMin,
                double tauMax,
                int numAveragingSteps )
    : BaseInfo ( TOL_Newton, maxIterNewton, solverInfo,
                 tauMin, tauMax, numAveragingSteps )
    , _wasInterrupted (false)
    , _previousCtrlCHandler ( ErrorSignal )
{}
//---------------------------------------------------------------------------

template <typename BaseInfo>
InterruptableIterativeInfo<BaseInfo>::InterruptableIterativeInfo
              ( RealType TOL_Newton,
                int maxIterNewton,
                RealType TOL_solver,
                int maxIterSolver,
                StoppingMode mode,
                ostream & out,
                double tauMin,
                double tauMax,
                int numAveragingSteps )
    : BaseInfo ( TOL_Newton, maxIterNewton,
                 *( new InterruptableIterativeInfo<SolverInfo<RealType> >
                           ( TOL_solver, maxIterSolver, mode, false, out ) ),
                 tauMin, tauMax, numAveragingSteps )
    , _wasInterrupted (false)
    , _previousCtrlCHandler ( ErrorSignal ) {
  this->_solverInfo.setDeleteFlag ( true );
}
//---------------------------------------------------------------------------

template <typename BaseInfo>
void InterruptableIterativeInfo<BaseInfo>::startIterations() {
  // install ctrl-C handler and save previous handler
  _previousCtrlCHandler = signal ( InterruptSignal, ctrlCHandler );
  _wasInterrupted = false;

  // call base class' startIterations().
  BaseInfo::startIterations();
}

//---------------------------------------------------------------------------

template <typename BaseInfo>
void InterruptableIterativeInfo<BaseInfo>::finishIterations() {
  if ( wasInterrupted() )
    this->_out << "Interrupted after " << this->getIterationCount()
               << " iterations." << endl << endl;

  BaseInfo::finishIterations();

  // reset ctrl-C handler that was in use
  // before iteration start
  if ( _previousCtrlCHandler != ErrorSignal )
    _previousCtrlCHandler = signal (InterruptSignal, _previousCtrlCHandler);
}

//---------------------------------------------------------------------------
template <typename BaseInfo>
bool InterruptableIterativeInfo<BaseInfo>::maxIterIsExceeded() const
             {  return (BaseInfo::maxIterIsReached() || wantsInterrupt());  }
//---------------------------------------------------------------------------
template <typename BaseInfo>
bool InterruptableIterativeInfo<BaseInfo>::maxIterIsReached() const
             {  return (BaseInfo::maxIterIsReached() || wantsInterrupt());  }
//---------------------------------------------------------------------------
template <typename BaseInfo>
bool InterruptableIterativeInfo<BaseInfo>::wasInterrupted() const
                                                {  return _wasInterrupted;  }
//---------------------------------------------------------------------------

template <typename BaseInfo>
bool InterruptableIterativeInfo<BaseInfo>::wantsInterrupt() const {
  if (!getCtrlCState())
    return false;

  signal ( InterruptSignal, DefaultHandler );
  cerr << endl << endl << aol::color::error
       << "Do you really want to interrupt the iteration "
          "(y/n, Ctrl-C to kill program)? ";
  cerr << aol::color::reset;
  string yes_no;
  cin >> yes_no;
  signal ( InterruptSignal, ctrlCHandler );

  if (yes_no == "y" || yes_no == "yes")
  {
    return (_wasInterrupted = true);
  }
  else
  {
    resetCtrlCState();
    return (_wasInterrupted = false);
  }
}

//---------------------------------------------------------------------------

} // end of namespace aol.

#endif

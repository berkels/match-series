#ifndef __PROGRESSBAR_H
#define __PROGRESSBAR_H

#include <aol.h>
#include <qmException.h>

namespace aol {

/** Class that shows a busy bar in ascii text
 *  \author Preusser
 */
template <int _screenWidth = 79, int _fullLength = 10, char _leftBrace = '[', char _rightBrace = ']',
char _empty = '.', char _full = '='>
class BusyBar {
protected:
  char _fullChar[_fullLength+1];
  char _emptyChar[_screenWidth];

  char *_frontCharPtr;
  char *_endCharPtr;

  int  _value, _delta, _dir, _barLength, _pos, _incr;

  const char *_text;
  std::ostream &_out;

  /** Create the string representation of the progress bar
   */
  void createBar ( int p ) {
    _frontCharPtr = _emptyChar + _screenWidth - 1 - p;
    _endCharPtr = _emptyChar + _fullLength + ( _screenWidth - _barLength ) + p ;
  }

public:
  /** Initialize a busy bar. The text is drawn in front of the bar. If
   *  no automatic output to out is desired, only thre prefix increment operator++
   *  should be used to update the bar.
   */
  BusyBar ( const char *text = NULL, std::ostream &out = std::cerr ) : _text ( text ), _out ( out ) {
    memset ( _fullChar, _full, _fullLength );
    _fullChar[_fullLength] = '\0';
    memset ( _emptyChar, _empty, _screenWidth - 1 );
    _emptyChar[_screenWidth-1] = '\0';

    _barLength = _screenWidth - 2;
    if ( text != NULL ) _barLength -= static_cast<int>( strlen ( _text ) ) + 1;
    start();
  }

  /** Start a busy bar which shifts by incr for each delta'th call of operator++
   *  At each call of operator ++ the internal counter is updated by incr
   *  every delta steps of the internal counter an updated output is made.
   */
  void start ( const int delta = 5, const int incr = 1 ) {
    _value = 0;
    _delta = delta;
    _pos = 0;
    _dir = 1;
    _incr = incr;
    createBar ( 0 );
  }

  /** Put all parts of the string representation together and show it.
   *  To avoid flickering of the cursor, we print into a stringstream first.
   */
  void display ( std::ostream &out ) const {
    std::ostringstream oss;
    oss << '\r';
    if ( _text ) oss << _text;
    oss << ' ' << _leftBrace << _frontCharPtr << _fullChar << _endCharPtr << _rightBrace << std::flush;
    out << oss.str();
  }

  /** postfix-increment: increment the internal counter and print progress bar afterwards
   */
  void operator++ ( int ) {
    if ( this->operator++() ) display ( _out );
  }

  /** prefix-increment: increment the internal counter but do not print the progress bar
   *  Returns true if an update was made
   */
  bool operator++() {
    _value += _incr;
    if ( _value % _delta == 0 ) {
      _pos += _dir;
      createBar ( _pos );
      if ( _pos == 0 || _pos == _barLength - _fullLength ) _dir *= -1;
      return true;
    }
    return false;
  }
};



/** Class that creates a progress bar in ascii text.
 *  \attention ProgressBars are not thread-safe (a "critical" section in the increment operator would make things extremely slow).
 *  \author Preusser
 */
template <bool _showPercentage = true, bool _showBar = true, typename TYPE = int, int _screenWidth = 79,
char _leftBrace = '[', char _rightBrace = ']', char _empty = '.', char _full = '='>
class ProgressBar {
protected:
  char _prefix;
  char _postfix;
  char _percent[10];
  char _fullChar[_screenWidth];
  char _emptyChar[_screenWidth];
  char _backSpace[_screenWidth];

  char *_fullCharPtr;
  char *_emptyCharPtr;
  char *_bsPointer;

  TYPE _value, _maxValue, _delta, _incr;

  const char *_text;
  std::ostream &_out;

  int    _barLength, _sw1, _swBl, _lastUpdatedQ;

  /** Create the string that shows the percentage
   */
  void createPercentage ( int p ) {
    if ( _showPercentage ) {
      sprintf ( _percent, " %3d%% ", p );
    }
  }

  /** Create the string representation of the progress bar
   */
  void createBar ( double p ) {
    if ( p > 1 ) // do not show more than 100 % in the string representation (else strings too short)
      p = 1;

    if ( _showBar ) {
      int q = static_cast<int> ( p * _barLength );
      const int v = _sw1 - q;
      _fullCharPtr = _fullChar + v;

      const int w = _swBl + q ;
      _emptyCharPtr = _emptyChar + w;
      _bsPointer = _backSpace + w - 1;
    }
  }

  void initConstants() {
    _barLength = _screenWidth;
    if ( _text != NULL ) _barLength -= static_cast<int> ( strlen ( _text ) );
    if ( _showPercentage ) _barLength -= 6;
    if ( _showBar ) {
      _prefix = _leftBrace;
      _postfix = _rightBrace;
      _barLength -= 2;
    } else {
      _prefix = '\0';
      _postfix = '\0';
    }

    _sw1 = _screenWidth - 1;
    _swBl = _screenWidth - 1 - _barLength;
  }
public:
  /** Initialize a progress bar. Display of the bar and the percentage can be turned
   *  off with the corresponding variables. The text is drawn in front of the bar. If
   *  no automatic output to out is desired, only thre prefix increment operator++
   *  should be used to update the bar.
   */
  ProgressBar ( const char *text = NULL, std::ostream &out = std::cerr ) : _text ( text ), _out ( out ) {
    if ( !_showPercentage && !_showBar ) throw aol::Exception ( "ProgressBar must show bar or percentage", __FILE__, __LINE__ );
    memset ( _backSpace, '\b', _screenWidth - 1 );
    _backSpace[_screenWidth-1] = '\0';
    memset ( _fullChar, _full, _screenWidth - 1 );
    _fullChar[_screenWidth-1] = '\0';
    memset ( _emptyChar, _empty, _screenWidth - 1 );
    _emptyChar[_screenWidth-1] = '\0';

    if ( !_showPercentage )
      cerr << aol::color::red << "aol::ProgressBar with first template parameter _showPercentage == false seems to be buggy" << aol::color::reset << endl;
    // remove this output if it has been fixed.

    initConstants();
  }
  
  virtual ~ProgressBar ( ) { }

  virtual void setText ( const char *text ) {
    _text = text;
    initConstants();
  }

  /** Start a progress with maxValue and an update for each delta Percent
   *  At each call of operator ++ the internal counter is updated by incr
   *  Every delta steps of the internal counter an updated output is made.
   */
  virtual void start ( const TYPE maxValue, const int delta = 5, const TYPE incr = 1 ) {
    _value = 0;
    _maxValue = maxValue;
    _delta = delta;
    _incr = incr;
    _lastUpdatedQ = -1;

    createPercentage ( 0 );
    createBar ( 0 );
  }
  void finish () const {
    _out << endl;
  }

  /** Display a progress of p percent
   */
  void set ( const int p = 100 ) {
    createPercentage ( p );
    createBar ( static_cast<double> ( p / 100.0 ) );
  }

  /** postfix-increment: increment the internal counter and print progress bar afterwards
   */
  virtual void operator++ ( int ) {
    if ( this->operator++() ) display ( _out );
  }

  /** prefix-increment: increment the internal counter but do not print the progress bar
   *  Returns true if an update was made
   */
  virtual bool operator++() {
    _value += _incr;
    const double p = ( _value * 1.0 / _maxValue );
    const int q = static_cast<int> ( p * 100 );
    if ( q != _lastUpdatedQ && q % _delta == 0 ) {
      createPercentage ( q );
      createBar ( p );
      _lastUpdatedQ = q;
      return true;
    }
    return false;
  }

  /** Put all parts of the string representation together and show it.
   *  To avoid flickering of the cursor, we print into a stringstream first.
   */
  virtual void display ( std::ostream &out ) const {
    std::ostringstream oss;
    oss << '\r' << _text << _percent << _prefix << _fullCharPtr << _emptyCharPtr << _postfix << _bsPointer << std::flush;
    out << oss.str();
  }
  
  const char* getPercentStr ( ) const {
    return _percent;
  }
};

/** A simplified progress bar estimating the remaining time (assuming a linear scale) in CPU time.
 *  This should work correctly when parallelizing (as wall clock time is used for the estimates), pausing , however, is not treated correctly.
 *  \author Schwen (MEVIS)
 */
class EstimatingProgressBar : public ProgressBar<> {
protected:
  double _startedAt;

  //! can be overloaded on derived classes to allow estimating nonlinear run time
  virtual double complexityEstimate ( const int Val ) const {
    return ( static_cast<double> ( Val ) );
  }

public:
  explicit EstimatingProgressBar ( const char* const text ) : ProgressBar<> ( text ), _startedAt ( 0.0 ) {
  }

  EstimatingProgressBar ( const char* const text, const int MaxValue ) : ProgressBar<> ( text ), _startedAt ( 0.0 ) {
    start ( MaxValue );
  }

  virtual ~EstimatingProgressBar ( ) {
    // no need to do anything, but destructor should be virtual since we have a virtual method
  }

  void start ( const int maxValue, const int delta = 1, const int /*incr*/ = 1 ) {
    aol::ProgressBar<>::start ( maxValue, delta );
    time_t nowTime;
    unsigned short nowMillisec;
    getWallClockTime( nowTime, nowMillisec );
    _startedAt = nowTime + 0.001 * nowMillisec;
  }

  void operator++ ( int ) {
    if ( operator++() ) display ( _out );
    const double p = ( _value * 1.0 / _maxValue );
    const int q = static_cast<int> ( p * 100 );
    if ( ( q != this->_lastUpdatedQ ) && ( q % _delta == 0 ) ) {
      time_t nowTime;
      unsigned short nowMillisec;
      getWallClockTime( nowTime, nowMillisec );

      const double
        dp = ( complexityEstimate ( _value ) / complexityEstimate ( _maxValue ) ),
        nowTimeD = nowTime + 0.001 * nowMillisec,
        elap = nowTimeD - _startedAt,
        estm = elap * ( 1 - dp ) / dp;

      const time_t estimFinished = static_cast<time_t> ( nowTimeD + estm );
      char etaString[512];
      strftime ( etaString, 512, "%a %F %T", localtime ( &estimFinished ) );

      _out << _emptyCharPtr << _postfix << " ETA " << etaString << ", " << aol::intFormat ( static_cast<int> ( estm ) ) << " s remaining.";
      this->_lastUpdatedQ = q;
    }
  }

  bool operator++() {
    _value += _incr;
    const double p = ( _value * 1.0 / _maxValue );
    const int q = static_cast<int> ( p * 100 );
    if ( ( q != this->_lastUpdatedQ ) && ( q % _delta == 0 ) ) {
      createPercentage ( q );
      createBar ( p );
      return true;
    }
    return false;
  }


  void finish ( ) const {
#ifdef VERBOSE
    const double elap = aol::getRuntimeSoFar() - _startedAt;
    _out << endl << "Took " << aol::longIntFormat ( static_cast<int> ( elap ) ) << " seconds.";
#endif
    aol::ProgressBar<>::finish();
  }
  
  const char* getTimeRemainingStr ( ) {
    time_t nowTime;
    unsigned short nowMillisec;
    getWallClockTime( nowTime, nowMillisec );
    const double
      dp = ( complexityEstimate ( _value ) / complexityEstimate ( _maxValue ) ),
      nowTimeD = nowTime + 0.001 * nowMillisec,
      elap = nowTimeD - _startedAt,
      estm = elap * ( 1 - dp ) / dp;
    const int
      days = estm / 86400,
      hours = ( estm - days * 86400 ) / 3600,
      minutes = ( estm - days * 86400 - hours * 3600 ) / 60,
      seconds = estm - days * 86400 - hours * 3600 - minutes * 60;
    std::stringstream ss;
    ss << "ETA ";
    if ( estm > 31536000 ) ss << "> 1 year";
    else {
      if ( days < 0 ) ss << days << "d ";
      if ( hours < 10 ) ss << "0";
      ss << hours << ":";
      if ( minutes < 10 ) ss << "0";
      ss << minutes << ":";
      if ( seconds < 10 ) ss << "0";
      ss << seconds;
    }
    return ss.str ( ).c_str ( );
  }
};


/** Output the busy bar
 */
template <int _screenWidth, int _fullLength, char _leftBrace, char _rightBrace,
char _empty, char _full>
inline std::ostream &operator<< ( std::ostream &out,
                                  const BusyBar < _screenWidth, _fullLength, _leftBrace, _rightBrace,
                                  _empty, _full > &bb ) {
  bb.display ( out );
  return out;
}

/** Output the progress bar
 */
template <bool _showPercentage, bool _showBar, typename TYPE, int _screenWidth,
char _leftBrace, char _rightBrace, char _empty, char _full>
inline std::ostream &operator<< ( std::ostream &out,
                                  const ProgressBar < _showPercentage, _showBar, TYPE, _screenWidth,
                                  _leftBrace, _rightBrace, _empty, _full > &pb ) {
  pb.display ( out );
  return out;
}

}

#endif

#ifndef __PARAMETERPARSER_H
#define __PARAMETERPARSER_H

#include <aol.h>
#include <vec.h>
#include <multiVector.h>

namespace aol {

class Variable {
public:
  enum { VAR_DOUBLE = 1,
         VAR_INT    = 2,
         VAR_STRING = 3 };

  Variable ( const char *VarStr ) {
    varStr = VarStr;
  }

  Variable() {
    varStr = "";
  }

  int type() const {
    return type ( getVarStr() );
  }

  // This comparison operator returning non-bool is not used.
  // int operator< ( const Variable &Var ) const {
  //  return strcmp ( Var.varStr, varStr );
  // }

  bool operator== ( const Variable &Var ) const {
    return !strcmp ( Var.getVarStr(), getVarStr() );
  }

  Variable &operator= ( const Variable &Var ) {
    varStr = Var.varStr;
    return *this;
  }

  void dump ( ostream &Off = cerr ) const {
    Off << varStr;
  }

  double getDouble() const {
    return atof ( getVarStr() );
  }

  int getInt() const {
    return atoi ( getVarStr() );
  }

  const char *getVarStr() const {
    return varStr.c_str();
  }

  void setVarStr ( const char *VarStr ) {
    varStr = VarStr;
  }

private:

// turning warnings off right before calling the critical function does only work with gcc > 4.5
WARNING_OFF ( unused-result )
  int type ( const char *VarStr ) const {

    char *err;
    //long int ivalue = strtol(VarStr, &err, 0);
// WARNING_OFF ( unused-result )
    strtol ( VarStr, &err, 0 );
// WARNING_ON ( unused-result )
    if ( *err == '\0' ) return VAR_INT;

    //double dvalue = strtod(VarStr, &err);
// WARNING_OFF ( unused-result )
    strtod ( VarStr, &err );
// WARNING_ON ( unused-result )
    if ( *err == '\0' ) return VAR_DOUBLE;

    return VAR_STRING;

  }
WARNING_ON ( unused-result )

  string varStr;
};

class VariableField {
public:
  VariableField ( const char *Name ) :
      numDim ( 0 ),
      dimSizes ( NULL ) {
    name =  Name;
  }

  VariableField ( int Num ) : dimSizes ( NULL ) {
    setNumDim ( Num );
  }

  VariableField ( const VariableField &Other )
    : dimSizes ( NULL ),
      name ( Other.name ),
      field ( Other.field.size() ) {
    setNumDim ( Other.numDim );

    for ( int i = 0; i < numDim; ++i )
      dimSizes[i] = Other.dimSizes[i];

    for ( unsigned int i = 0; i < Other.field.size(); ++i )
      field[i] = new Variable ( *(Other.field[i]) );
  }

  ~VariableField() {
    if ( dimSizes )
      delete[] dimSizes;
    vector<Variable*>::const_iterator it;
    for ( it = field.begin(); it != field.end(); ++it ) {
      delete *it;
    }
  }

  void setNumDim ( int Num ) {
    numDim = Num;
    if ( dimSizes )
      delete[] dimSizes;
    dimSizes = new int[ Num ];
  }

  int getNumDim() const {
    return numDim;
  }

  void setDimSize ( int Num, int Size ) {
    if ( Num < numDim && Num >= 0 ) {
      dimSizes[ Num ] = Size;
    } else {
      throw Exception ( "ERROR in VariableField::setDimSize: Num out of range.\n", __FILE__, __LINE__ );
    }
  }

  int getDimSize ( int Num ) const {
    if ( Num < numDim && Num >= 0 ) {
      return dimSizes[ Num ];
    } else {
      throw Exception ( "ERROR in VariableField::getDimSize: Num out of range.\n", __FILE__, __LINE__ );
      return -1;
    }
  }

  void append ( const char *Var ) {
    /*cerr << "VarField::append( " << Var << " ) \n" ;*/
    field.push_back ( new Variable ( Var ) );
  }

  void dump ( ostream &Off = cerr ) const {
    Off << "VarField: " << name << endl;
    Off << "numDim = " << numDim << endl;
    vector<Variable*>::const_iterator it;
    for ( it = field.begin(); it != field.end(); ++it ) {
      ( *it )->dump ( Off );
      Off << ", ";
    }
    Off << endl;
    if ( numDim > 0 ) {
      for ( int i = 0; i < numDim; i++ ) {
        Off << "dimSizes[ " << i << " ]= " << dimSizes[ i ] << endl;
      }
    }
  }

  //! Write contents to a stream in a format parsable by aol::ParameterParser
  //! \author Berkels
  void write ( ostream &Off = cerr ) const;

  int isSingleField() const {
    return ( numDim == 0 );
  }

  void getVariable ( int I, Variable &Var ) const {
    if ( I < static_cast<int> ( field.size() ) ) {
      Var = *field[ I ];
    } else {
      cerr << "index out of bounds\n";
    }
  }

  void getVariable ( int I1, int I2, Variable &Var ) const {
    if ( I1 < dimSizes[ 0 ] && I2 < dimSizes[ 1 ] ) {
      Var = *field[ I1 * dimSizes[ 1 ] + I2 ];
    } else {
      cerr << "index out of bounds\n";
    }
  }

  void getVariable ( Variable &Var ) const {
    if ( field.size() == 1 ) {
      Var = *field[ 0 ];
    } else {
      cerr << "ERROR in VariableField::getVariable: size = " << field.size() << endl;
    }
  }

  vector<Variable*> &getFieldReference ( ) {
    return field;
  }

  const char *getName() const {
    return name.c_str();
  }

protected:
  int numDim;
  int *dimSizes;
  string name;
  vector<Variable*> field;
};

/**
 * \brief A class for reading parameters from configuration files.
 *
 * A parameter file may contain the following type of lines:
 * - parameter line, i.e. PARAMETERNAME VALUE
 * - comment lines, i.e. lines starting with #
 * - white space lines, i.e. lines containing only white spaces
 * - blank lines, i.e. lines containing only EOL
 *
 * If the file contains any other type of lines, expect any kind of odd behavior.
 * \author Droske
 *
 * Currently maintained by Berkels.
 *
 * \todo Rewrite all the functions that look for a variable to use findFirstVariableField.
 * \ingroup Utilities
 */
class ParameterParser {
public:
  //! Default constructor. Since this doesn't parse any parameter file,
  //! you can't do much with an instance created by this constructor.
  ParameterParser ( ) : _echo ( false ), _outStream ( std::cout ) { }

  /*! Instantiate a parser for a certain file. */
  explicit ParameterParser ( const std::string & ParFilename ) : _echo ( false ), _outStream ( std::cout ) {
    initialize ( ParFilename );
  }

  /*! Instantiate a parser for the file named argv[1] if argc == 2 or named DefaultParFilename if argc == 1. */
  ParameterParser ( int argc, char **argv, const char *DefaultParFilename );

  //! creates a ParameterParser from fully constructed other ParameterParser
  //! does therefore no copying of ifstream in
  ParameterParser ( const ParameterParser & other );

  ~ParameterParser() {
    vector<VariableField*>::iterator it;
    for ( it = varFields.begin(); it != varFields.end(); ++it ) {
      delete *it;
    }
  }


  int    getDimSize ( const char *VarName, int I = 0 ) const;
  int    getNumDim ( const char *VarName ) const;

  double getDouble ( const char *VarName ) const;
  double getDouble ( const char *VarName, int I ) const;
  double getDouble ( const char *VarName, int I1, int I2 ) const;
  double getDoubleOrDefault ( const char *VarName, double Default ) const;

  template <typename RealType>
  RealType getReal ( const char *VarName ) const {
    return static_cast<RealType>( this->getDouble ( VarName ) );
  }

  template <typename RealType>
  RealType getRealOrDefault ( const char *VarName, double Default ) const {
    return static_cast<RealType>( this->getDoubleOrDefault( VarName, Default ) );
  }

  int    getInt ( const char *VarName, int I ) const;
  int    getInt ( const char *VarName, int I1, int I2 ) const;
  int    getInt ( const char *VarName ) const;
  int    getIntOrDefault ( const char *VarName, int Default ) const;

  void   getString ( const char *VarName, char *DestStr ) const;
  void   getString ( const char *VarName, char *DestStr, int I ) const;
  string getString ( const char *VarName ) const;
  string getString ( const char *VarName, const int I ) const;
  string getStringOrDefault ( const char *VarName, string Default ) const;

  string getStringExpandTilde ( const char *VarName ) const;

  //! check whether parameter file has a parameter
  bool   hasVariable ( const char* VarName ) const;

  //! check whether parameter file has a variable, print
  //! error message to cerr if not
  bool   checkVariable ( const char * VarName ) const;

  //! Returns true if the parameter file has the variable \arg VarName and its int value is 1.
  //! Otherwise, false is returned.
  bool   checkAndGetBool ( const char *VarName ) const;

  //! Returns true/false if the parameter file has the variable \arg VarName and its int value is 1/0.
  //! Throws an exception if either the variable doesn't exist or has a value different from 0 and 1.
  bool   getBool ( const char *VarName ) const;
  
  //! Returns true/false if the parameter file has the variable \arg VarName and its int value is 1/0
  //! Otherwise, returns \arg Default
  bool getBoolOrDefault ( const char *VarName, bool Default ) const;

  void   dump ( ostream &Off = cerr ) const;
  void   dumpToFile ( const char *FileName, const char *Directory = NULL ) const;

  void   setEchoMode ( const bool echo = true ) {
    _echo = echo;
  }

  string getParameterFileName () const {
    return _parameterFileName;
  }

  template<typename RealType>
  void getRealVec ( const char *VarName, aol::Vector<RealType> &vec ) const {
    const int size = this->getDimSize ( VarName );
    vec.resize ( size ); vec.setZero();
    for ( int i = 0; i < size; ++i ) {
      vec[i] = static_cast<RealType>( this->getDouble ( VarName, i ) );
    }
  }

  template<int NumComponents, typename RealType>
  aol::Vec<NumComponents, RealType> getRealVec ( const char *VarName ) const {
    if ( this->getDimSize ( VarName ) != NumComponents )
      throw Exception ( "Parameter entry has wrong size.\n", __FILE__, __LINE__ );

    aol::Vec<NumComponents, RealType> vec;
    for ( int i = 0; i < NumComponents; ++i ) {
      vec[i] = static_cast<RealType>( this->getDouble ( VarName, i ) );
    }
    return vec;
  }

  //! Assumes that the field is given as <tt>{ {a^1_1 ... a^1_n} ... {a^m_1 ... a^m_n} }</tt>
  template<typename RealType>
  void getRealMultiVec ( const char *VarName, aol::MultiVector<RealType> &MVec ) const {
    if ( this->getNumDim ( VarName ) != 2 )
      throw Exception ( "Field has an improper format.\n", __FILE__, __LINE__ );

    const int numVec = this->getDimSize ( VarName, 0 );
    const int vecSize = this->getDimSize ( VarName, 1 );
    MVec.reallocate ( numVec, vecSize );
    for ( int i = 0; i < numVec; ++i )
      for ( int j = 0; j < vecSize; ++j )
        MVec[i][j] = static_cast<RealType>( this->getDouble ( VarName, i, j ) );
  }

  template<typename RealType>
  void getRealMultiVecTransposed ( const char *VarName, aol::MultiVector<RealType> &MVec ) const {
    if ( this->getNumDim ( VarName ) != 2 )
      throw Exception ( "Field has an improper format.\n", __FILE__, __LINE__ );

    const int numVec = this->getDimSize ( VarName, 0 );
    const int vecSize = this->getDimSize ( VarName, 1 );
    MVec.reallocate ( vecSize, numVec );
    for ( int i = 0; i < numVec; ++i )
      for ( int j = 0; j < vecSize; ++j )
        MVec[j][i] = static_cast<RealType>( this->getDouble ( VarName, i, j ) );
  }

  void   getIntVec ( const char *VarName, aol::Vector<int> &vec ) const;

  template<int NumComponents>
  aol::Vec<NumComponents, int> getIntVec ( const char *VarName ) const {
    if ( this->getDimSize ( VarName ) != NumComponents )
      throw Exception ( "Parameter entry has wrong size.\n", __FILE__, __LINE__ );

    aol::Vec<NumComponents, int> vec;
    for ( int i = 0; i < NumComponents; ++i )
      vec[i] = this->getInt ( VarName, i );
    return vec;
  }

  void changeVariableValue ( const char *VarName, const char *NewValue ) {
    VariableField* var = findFirstVariableField ( VarName );

    if ( var->isSingleField() == false )
      throw DimensionMismatchException ( "VariableField::changeVariableValue only implemented for single fields.\n", __FILE__, __LINE__ );

    var->getFieldReference()[0]->setVarStr ( NewValue );
  }

  void changeVariableValue ( const char *VarName, const int NewValue ) {
    changeVariableValue ( VarName, aol::strprintf ( "%d", NewValue ).c_str() );
  }
  
  template<typename RealType>
  void changeVariableValue ( const char *VarName, const RealType NewValue ) {
    changeVariableValue ( VarName, aol::strprintf ( "%f", NewValue ).c_str() );
  }

  void addVariable ( const char *VarName, const char *Value ) {
    if ( hasVariable ( VarName ) )
      throw Exception ( aol::strprintf ( "Variable \"%s\" already exists.\n", VarName ).c_str(), __FILE__, __LINE__ );

    VariableField *varField = new VariableField ( VarName );
    varField->setNumDim ( 0 );
    varField->append ( Value );
    varFields.push_back ( varField );
  }

  void addVariable ( const char *VarName, const int NewValue ) {
    addVariable ( VarName, aol::strprintf ( "%d", NewValue ).c_str() );
  }

protected:
  VariableField* findFirstVariableField ( const char *VarName ) const;
  void initialize ( const std::string & ParFilename );
  void parse();
  void readField ( VariableField &VarField );

  int ws ( int c ) const {
    return ( ( c == ' ' || c <= 0x17 ) && ( c != '\n' ) );
  }

  void ignoreWS() {
    while ( !in.eof() && ws ( in.peek() ) ) {
      /*cerr << "ignoring << " << in.get( ) << endl;*/
      in.ignore();
    }
  }

  void ignoreLine() {
    while ( !in.eof() && in.peek() != '\n' ) {
      in.ignore();
    }
    if ( !in.eof() && in.peek() == '\n' )
      in.ignore();
  }

  vector<VariableField*> varFields;

  ifstream in;

  bool     _echo;
  ostream  &_outStream;
  string   _outString;
  string   _parameterFileName;
};

void addCounterToSaveDirectory ( ParameterParser&, const string );

/**
 * \author Berkels
 */
void createTemplateFileNameList ( const aol::ParameterParser &Parser, std::vector<std::string> &FileNames );

}

#endif

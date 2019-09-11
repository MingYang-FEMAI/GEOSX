/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef GEOSX_SRC_CORECOMPONENTS_FILEIO_LAS_LASFILE_
#define GEOSX_SRC_CORECOMPONENTS_FILEIO_LAS_LASFILE_

#include "common/DataTypes.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

class LASVersionInformationSection;

/*!
 * @brief This class contains a parsed line of a LAS file
 * @details A LAS line have always this structure
 * MNEM.UNIT    DATA  :  DESCRIPTION
 * With '.' as a first delimiter, the first space as the second
 * delimiter and ':' as the third delimiter
 */
class LASLine
{
  public:
    LASLine() = default;

    LASLine( string const& line )
    {
      ParseLine( line );
    }

    string const & GetKeyword() const
    {
      return m_keywordname;
    }

    string const & GetData() const
    {
      return m_data;
    }

    bool GetDataAsBool() const
    {
      GEOS_ASSERT( m_data == "YES" || m_data == "NO");
      return ( m_data == "YES" ? true : false );
    }

    localIndex GetDataAsLocalIndex() const
    {
      return stringutilities::fromString< localIndex >( m_data );
    }

    real64 GetDataAsReal64() const
    {
      return stringutilities::fromString< real64 >( m_data );
    }

    string const & GetDescription() const
    {
      return m_description;
    }

    string const & GetUnit() const
    {
      return m_unit;
    }

    bool HasUnit() const
    {
      return !m_unit.empty();
    }
  private:
    void ParseLine( string const & line )
    {
      // First get the keyword and the rest of the line
      string_array keywordAndRest = stringutilities::Split( line, "." );
      stringutilities::RemoveSpaces( keywordAndRest[0] );
      m_keywordname = keywordAndRest[0];

      // Second get the unit and the rest of the line
      string_array unitsAndRest = stringutilities::Split( keywordAndRest[1], " ");
      stringutilities::RemoveSpaces( unitsAndRest[0] );
      m_unit = unitsAndRest[0];

      // Third get the value and the rest of te line
      string_array valueAndRest = stringutilities::Split( unitsAndRest[1], " " );
      stringutilities::Trim( valueAndRest[0] );
      m_data = valueAndRest[0];

      // Finally, get the description
      stringutilities::Trim( valueAndRest[1] );
      m_description = valueAndRest[1];
    }
  private:
    /// Name of the keyword
    string m_keywordname;

    /// Unit (if applicable);
    string m_unit;

    /// Data
    string m_data;

    /// Description
    string m_description;
};

/*!
 * @brief Basis class for a LAS Information Section
 * @details a LAS Information Section contains an ensemble of LASLine,
 * referenced by their keywords
 */
class LASSection
{
  public:
    LASSection() = default;
    virtual ~LASSection()
    {
    }

    /*!
     * Check if all the mandatory keywords are registered
     */
    virtual void CheckKeywords()
    {
      for( string & keyword : m_mandatoryKeyword )
      {
        GEOS_ERROR_IF( m_lines.count( keyword ) == 0, "Mandatory keyword " << keyword << " not found in "
                                                      << GetName() );
      }
    }

    virtual std::streampos ParseSection( std::ifstream & file )
    {
      string curLine;
      std::streampos pos;
      while ( std::getline(file, curLine) )
      {
        stringutilities::TrimLeft( curLine );
        if( curLine[0] == '#' ) continue;
        if( curLine[0] == '~' )                // We reach a new section
        {
          break;
        }
        pos = file.tellg();
        ParseLine( curLine );
      }
      return pos;
    }
    
    virtual string const GetName() const = 0;

    static std::unique_ptr< LASSection > CreateSection( char const & name );

  protected:
    virtual void ParseLine( string const & line )
    {
      LASLine curLine( line );
      GEOS_ERROR_IF( m_lines.count( curLine.GetKeyword() ) != 0, "Keyword " << curLine.GetKeyword()
                     << " was already defined in "<< GetName() );
      m_lines[curLine.GetKeyword()] = curLine;
    }


  protected:
    /// Contains the mandatory keyword for the section
    array1d< string > m_mandatoryKeyword;

    /// Contains all the lines (A line = keyword + value(s))
    std::unordered_map< string, LASLine > m_lines;
};

/*!
 * @brief This section identify the LAS format ans whether the wrap format is used
 */
class LASVersionInformationSection : public LASSection
{
  public:
    LASVersionInformationSection() :
      LASSection()
    {
      m_mandatoryKeyword.push_back( "VERS" );
      m_mandatoryKeyword.push_back( "WRAP" );
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Version Information";
    }

};

/*!
 * @brief This section identify the well
 */
class LASWellInformationSection : public LASSection
{
  public:
    LASWellInformationSection() :
      LASSection()
    {
      m_mandatoryKeyword.push_back( "STRT" );
      m_mandatoryKeyword.push_back( "STOP" );
      m_mandatoryKeyword.push_back( "STEP" );
      m_mandatoryKeyword.push_back( "NULL" );
      m_mandatoryKeyword.push_back( "COMP" );
      m_mandatoryKeyword.push_back( "WELL" );
      m_mandatoryKeyword.push_back( "FLD" );
      m_mandatoryKeyword.push_back( "LOC" );
      m_mandatoryKeyword.push_back( "PROV" );
      m_mandatoryKeyword.push_back( "CNTY" );
      m_mandatoryKeyword.push_back( "STAT" );
      m_mandatoryKeyword.push_back( "CTRY" );
      m_mandatoryKeyword.push_back( "SRCV" );
      m_mandatoryKeyword.push_back( "DATE" );
      m_mandatoryKeyword.push_back( "UWI" );
      m_mandatoryKeyword.push_back( "API" );
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Well Information";
    }
};

/*!
 * @brief This section describes the curves and their units
 */
class LASCurveInformationSection : public LASSection
{
  public:
    LASCurveInformationSection() :
      LASSection()
    {
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Curve Information";
    }

    /*!
     * @brief Returns the number of curves
     * @details DEPTH, DEPT, TIME or INDEX is one the curve
     * A curve is a log
     */
    localIndex GetNumberOfCurves() const
    {
      return m_lines.size();
    }
  private:
    /*!
     * @brief For this specific section, one keyword must be
     * DEPTH, DEPT, TIME or INDEX"
     * */
    virtual void CheckKeywords() override
    {
      GEOS_ERROR_IF( m_lines.count( "DEPTH" ) + m_lines.count( "DEPT" ) +  m_lines.count( "TIME" ) +  m_lines.count( "INDEX ") != 1,
                     "Invalid " << GetName() << " section. It musts contains at least one those keyword \n"
                     << "DEPTH DEPT TIME INDEX");
    }
};

/*!
 * @brief This optional section describes optional parameters
 */
class LASParameterInformationSection : public LASSection
{
  public:
    LASParameterInformationSection() :
      LASSection()
    {
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Parameter Information";
    }
};

/*!
 * @brief This optional section describes some comments
 */
class LASOtherInformationSection : public LASSection
{
  public:
    LASOtherInformationSection() :
      LASSection()
    {
    }


    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Other Information";
    }
  private:
    virtual void ParseLine( string const & line ) override
    {
      m_comments += line;
    }

    /*!
     * @brief For this specific section there is no keyword
     * */
    virtual void CheckKeywords() override
    {
      GEOS_ERROR_IF( !m_lines.empty(),
                     "Invalid " << GetName() << " section. No keyword should be defined, only data.");
    }

  private:
    /// Content of the section
    string m_comments;

};

/*!
 * @brief This optional section contains all the log data
 */
class LASASCIILogDataSection : public LASSection
{
  public:
    LASASCIILogDataSection() :
      LASSection() ,
      m_logs(),
      m_curPos(0)
    {
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "ASCII Log Data";
    }
  private:
    /*!
     * @brief For this specific section there is no keyword
     * */
    virtual void CheckKeywords() override
    {
      GEOS_ERROR_IF( !m_lines.empty(),
                     "Invalid " << GetName() << " section. No keyword should be defined, only data.");
    }

    virtual std::streampos ParseSection( std::ifstream & file ) override
    {
      string curLine;
      std::streampos beginPos;
      // First pass to count the number of value
      localIndex nbValues = 0;
      while ( std::getline(file, curLine) )
      {
        string_array splitLine = stringutilities::Split( curLine, " " );
        if( curLine[0] == '#' ) continue;
        if( curLine[0] == '~' )                // We reach a new section
        {
          break;
        }
        nbValues++;
      }
      m_logs.resize( nbValues );

      file.seekg( beginPos );
      return LASSection::ParseSection( file );
    }

    virtual void ParseLine( string const & line ) override
    {
      string_array splitLine = stringutilities::Split( line, " " );
      for( string & linePart : splitLine )
      {
        m_logs[m_curPos++] = stringutilities::fromString< real64 >( linePart );
      }
    }
  private:
    /// Contains the well logs
    array1d< real64 > m_logs;

    /// Position on m_logs while reading the file
    localIndex m_curPos;
};

class LASFile
{
  public:
    LASFile( string const& fileName ):
      m_fileName( fileName )
    {
    }
    
    void Load()
    {
      std::ifstream file( m_fileName );
      GEOS_ERROR_IF( file.is_open(), "Can't open " << m_fileName );
      string curLine;
      while ( std::getline(file, curLine) )
      {
        stringutilities::TrimLeft( curLine );
        if( curLine[0] == '#' ) continue;  // Comment

        if( curLine[0] == '~' )            // Section
        {
          std::unique_ptr< LASSection > curLASSection = LASSection::CreateSection( curLine[1] );
          curLASSection->CheckKeywords();
          std::streampos curPos = curLASSection->ParseSection( file );
          file.seekg( curPos );
          m_lasSections[curLASSection->GetName()] = std::move( curLASSection );
        }
      }
    }
    
    void Save() const
    {
    }
  private:

    /// Name of the file to be loaded or save
    string m_fileName;

    /// Map containing all the LAS sections indexed by name
    std::unordered_map< string, std::unique_ptr< LASSection > >  m_lasSections;

};
}

#endif

/*
  Analysis Abstract class
*/

#include "CompactStar/Core/Analysis.hpp"

//==============================================================
//                        Analysis class
//==============================================================

/// @brief Constructor for the Analysis class.
/// Initializes the base Prog class with the label "Analysis".
CompactStar::Analysis::Analysis() 
: Prog("Analysis")
{ }

//--------------------------------------------------------------

/// @brief Virtual destructor for the Analysis class.
CompactStar::Analysis::~Analysis() { }

//--------------------------------------------------------------

/// @brief Sets the label for the analysis.
/// @param in_label A string label used to tag or identify this analysis.
void CompactStar::Analysis::SetLabel(const std::string& in_label) 
{
  label = in_label ;
}

//--------------------------------------------------------------

//==============================================================
/*
  Analysis Abstract class
*/

#include "CompactStar/Core/Analysis.hpp"

//==============================================================
//                        Analysis class
//==============================================================
// Constructor
CompactStar::Analysis::Analysis() 
: Prog("Analysis")
{ }
//--------------------------------------------------------------
/// Destructor
CompactStar::Analysis::~Analysis() { }
//--------------------------------------------------------------
// Sets the label
void CompactStar::Analysis::SetLabel(const std::string& in_label) 
{
  label = in_label ;
}
//--------------------------------------------------------------

//==============================================================

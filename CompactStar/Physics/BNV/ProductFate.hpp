#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace CompactStar::Physics::BNV
{

enum class TerminalProductFate { PromptEscape, SmThermalization, BoundInert, BoundInteracting };

struct ProductFateBranch
{
    std::string terminal_id;
    TerminalProductFate fate = TerminalProductFate::BoundInert;
    double weight = 0;
    std::string channel_id;
};

class ProductFateLedger final
{
  public:
    ProductFateLedger(std::string identity, std::vector<ProductFateBranch> branches)
      :identity_(std::move(identity)),branches_(std::move(branches))
    {
        if (identity_.empty() || branches_.empty()) throw std::runtime_error("missing product fate");
        std::set<std::pair<std::string,std::string>> ids;
        std::map<std::string,double> sums;
        for (const auto& b:branches_)
        {
            if (b.terminal_id.empty() || b.channel_id.empty() || !ids.insert({b.channel_id,b.terminal_id}).second)
                throw std::runtime_error("duplicate terminal product fate");
            if (!(b.weight>=0) || !std::isfinite(b.weight)) throw std::runtime_error("invalid product fate weight");
            sums[b.channel_id]+=b.weight;
        }
        for(const auto& sum:sums){const double budget=32*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(sum.second));
          if(std::abs(sum.second-1)>budget)throw std::runtime_error("product fate weights do not sum to one");}
    }
    const std::string& Identity() const { return identity_; }
    const std::vector<ProductFateBranch>& Branches() const { return branches_; }
    void RequireChannel(const std::string& channel)const
    {if(channel.empty()||std::none_of(branches_.begin(),branches_.end(),[&](const auto& b){return b.channel_id==channel;}))throw std::runtime_error("unowned product-fate channel");}
  private:
    std::string identity_;
    std::vector<ProductFateBranch> branches_;
};

} // namespace CompactStar::Physics::BNV

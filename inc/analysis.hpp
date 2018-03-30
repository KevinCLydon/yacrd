/*
Copyright (c) 2018 Pierre Marijon <pierre.marijon@inria.fr>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

/* standard include */
#include <memory>
#include <string>

/* project include */
#include "utils.hpp"

namespace yacrd {
namespace analysis {

void find_chimera(const std::string& paf_filename, std::uint64_t coverage_min, std::unordered_set<std::string>* remove_reads);

void add_gap(std::vector<yacrd::utils::interval>& middle, std::vector<yacrd::utils::interval>& extremity, yacrd::utils::interval& gap, std::uint64_t readlen);

} // namespace analysis
} // namespace yacrd

#endif // ANALYSIS_HPP

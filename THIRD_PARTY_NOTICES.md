# Third-party notices

The browser distribution contains transformed ANARCI model data and an independent Rust implementation of behavior described by ANARCI and HMMER. It also contains transformed germline data originating from IMGT. No native ANARCI, HMMER, MUSCLE, or Python code executes at runtime.

## Rust SIMD support

The runtime uses `wide` 1.6.1, `safe_arch` 1.2.0, and `bytemuck` 1.25.2. These crates are copyright Daniel "Lokathor" Gee and are available under Zlib OR Apache-2.0 OR MIT; this distribution uses them under Apache-2.0. Their source repositories are <https://github.com/Lokathor/wide>, <https://github.com/Lokathor/safe_arch>, and <https://github.com/Lokathor/bytemuck>.

## ANARCI

Pipeline and behavioral source: the official ANARCI repository, <https://github.com/oxpig/ANARCI>, commit `79f6c575056dedef86cb8f405ebb039197923eec`. The embedded profile and germline inputs were regenerated with that source pipeline from the IMGT/GENE-DB 3.1.43 live reference directory retrieved on 2026-09-23.

Copyright 2019 Charlotte Deane, James Dunbar, Alexsandr Kovaltsuk, Claire Marks

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Changes: `ALL.hmm` emissions, transitions, local-entry scores, and calibration values were converted into a browser-specific binary form; scores are quantized at 1/32768 natural-log units in signed 24-bit fields. The aligned V/J tables in `germlines.py` were converted to a five-bit packed alphabet. ANARCI's Python code is not bundled.

## HMMER

Reference source: HMMER 3.4, tag `hmmer-3.4`.

HMMER - Biological sequence analysis with profile hidden Markov models

Copyright (C) 1992-2023 Sean R. Eddy  
Copyright (C) 2015-2023 President and Fellows of Harvard College  
Copyright (C) 2000-2023 Howard Hughes Medical Institute  
Copyright (C) 1995-2006 Washington University School of Medicine  
Copyright (C) 1992-1995 MRC Laboratory of Molecular Biology

HMMER source code is distributed as open source under the terms of the BSD three-clause license:

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of any copyright holder nor the names of contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This project does not compile or bundle the HMMER or Easel C libraries. The Plan7 local/multihit recurrences used by the Rust implementation were independently expressed with HMMER 3.4 as the behavioral reference. HMMER's patented SIMD implementation is not copied or used.

## IMGT data and numbering

The germline sequences underlying ANARCI's models and aligned germline table originate from IMGT®, the international ImMunoGeneTics information system®, <https://www.imgt.org>. IMGT data and metadata are provided under the [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).

Attribution: IMGT®, the international ImMunoGeneTics information system®, Institute of Human Genetics (Université de Montpellier and CNRS), <https://www.imgt.org>.

Changes: the data were retrieved from the IMGT/GENE-DB live reference directory through ANARCI's pinned build pipeline on 2026-09-23. The pipeline curated and aligned the records and built the HMMER profiles; this project then packed the aligned records and quantized model scores as described above. No IMGT tool executes at runtime. IMGT®, Université de Montpellier, and CNRS do not endorse this project. IMGT® is a registered mark of CNRS.

# Third-party notices

Project-specific additions and modifications are distributed under the GNU General Public License, version 3 or, at your option, any later version (`GPL-3.0-or-later`). See [LICENSE](LICENSE). Third-party code retains the notices and permissions below; this document does not replace those notices.

## NSGA-II implementation

The following files retain the original copyright and BSD 2-Clause terms of Aravind Seshadri's implementation. Some files have been modified for this project's trajectory objectives, dynamic decision variables, checkpoint handling and population diversity.

- `SR5NURBS/SR5NURBS/NSGA2.m`
- `SR5NURBS/SR5NURBS/genetic_operator.m`
- `SR5NURBS/SR5NURBS/initialize_variables.m`
- `SR5NURBS/SR5NURBS/non_domination_sort_mod.m`
- `SR5NURBS/SR5NURBS/replace_chromosome.m`
- `SR5NURBS/SR5NURBS/tournament_selection.m`

The original license text is preserved in each file and reproduced here:

```text
Copyright (c) 2009, Aravind Seshadri
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in
     the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

## NURBS knot-span compatibility function

`SR5NURBS/SR5NURBS/findspan.m` is an unmodified copy of the function shipped with the NURBS toolbox. It is retained for compatibility; the current main pipeline uses its own basis and curve evaluation functions.

```text
Copyright (C) 2010 Rafael Vazquez

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
```

The repository's [LICENSE](LICENSE) supplies the GPL v3 text. The original copyright and license header in `findspan.m` is preserved.

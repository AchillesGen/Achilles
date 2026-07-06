// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef PDFBASE_HH
#define PDFBASE_HH

namespace achilles {

class PDFBase {
  public:
    virtual ~PDFBase() = default;
    virtual double operator()(double, double) const = 0;
};

} // namespace achilles

#endif

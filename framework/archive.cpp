///////////////////////////////////////////////////////////////////////////////////////
/// \file archive.cpp
/// \brief Classes to make (de)serializing to/from streams convenient
///
/// $Date: 2023-04-17 18:15:07 +0200 (Mon, 17 Apr 2023) $
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "archive.h"

ArchiveInStream::ArchiveInStream(std::istream& strm)
	: in(strm) {
}

bool ArchiveInStream::save() const {
	return false;
}

void ArchiveInStream::transfer(char* s, std::streamsize n) {
	in.read(s, n);
}

ArchiveOutStream::ArchiveOutStream(std::ostream& strm)
	: out(strm) {
}

bool ArchiveOutStream::save() const {
	return true;
}

void ArchiveOutStream::transfer(char* s, std::streamsize n) {
	out.write(s, n);
}

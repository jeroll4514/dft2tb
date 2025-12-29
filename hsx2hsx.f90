! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program hsx2hsx

! Generates an HS file from an HSX file (more compact)

use hsx_m, only: read_hsx_file, hsx_t, write_hsx_file

type(hsx_t)  :: h

call read_hsx_file(h,"HSX") 

call write_hsx_file(h, "HSX_out")

end program hsx2hsx


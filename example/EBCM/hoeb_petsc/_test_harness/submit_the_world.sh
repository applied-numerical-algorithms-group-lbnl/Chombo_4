#!/bin/csh -f

/submit.all_regular_debug_2d.jobs  _merged_debug_3d_6_28_2023/submit.merged_debug_3d.jobs
_all_regular_debug_3d_6_28_2023/submit.all_regular_debug_3d.jobs  _unmerged_debug_2d_6_28_2023/submit.unmerged_debug_2d.jobs
_merged_debug_2d_6_28_2023/submit.merged_debug_2d.jobs		  _unmerged_debug_3d_6_28_2023/submit.unmerged_debug_3d.jobs

echo "cd _all_regular_debug_2d_6_28_2023; source ./submit.all_regular_debug_2d.jobssubmit.all_regular_2d.jobs >& submit.out; cd .."
cd       _all_regular_2d_6_2_2023; source ./submit.all_regular_2d.jobs >& submit.out; cd ..

echo "cd      _merged_2d_6_2_2023; source      ./submit.merged_2d.jobs >& submit.out; cd .."
cd            _merged_2d_6_2_2023; source      ./submit.merged_2d.jobs >& submit.out; cd ..

echo "cd    _unmerged_2d_6_2_2023; source    ./submit.unmerged_2d.jobs >& submit.out; cd .."
cd          _unmerged_2d_6_2_2023; source    ./submit.unmerged_2d.jobs >& submit.out; cd ..

echo "cd _all_regular_3d_6_2_2023; source ./submit.all_regular_3d.jobs >& submit.out; cd .."
cd       _all_regular_3d_6_2_2023; source ./submit.all_regular_3d.jobs >& submit.out; cd ..

echo "cd      _merged_3d_6_2_2023; source      ./submit.merged_3d.jobs >& submit.out; cd .."
cd            _merged_3d_6_2_2023; source      ./submit.merged_3d.jobs >& submit.out; cd ..

echo "cd    _unmerged_3d_6_2_2023; source    ./submit.unmerged_3d.jobs >& submit.out; cd .."
cd          _unmerged_3d_6_2_2023; source    ./submit.unmerged_3d.jobs >& submit.out; cd ..

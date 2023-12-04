## error

Fatal error:
The group cutoff scheme has been removed since GROMACS 2020. Please use the
Verlet cutoff scheme.

## Solution

init.mdpの以下の部分を変更 

### 変更前
```
cutoff_scheme            = group
```
### 変更後
```
cutoff_scheme            = verlet  
```
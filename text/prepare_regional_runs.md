# ToDo’s until Tuesday 4th of August

- [ ]  Log into Bremen cluster
- [ ]  Prepare equilibrium runs for constant and random climate for Rofental
- [ ]  Create (empty) structure/layout of thesis in Latex

## Cluster Bremen:

- I don’t have access to the email with the instructions from Timo anymore (UIBK employee accounts).
- Asked Fabi to forward me the email…

## Rofental runs

The following glaciers make problems:

- Random climate with no temperature bias

  - RGI60-11.00737

    ```
    /Users/oberrauch/oggm-fork/oggm/core/vascaling.py:1737: RuntimeWarning: overflow encountered in double_scalars
      self.tau_a = self.tau_l * self.area_m2 / self.length_m ** 2
    /Users/oberrauch/oggm-fork/oggm/core/vascaling.py:1782: RuntimeWarning: divide by zero encountered in double_scalars
      - self.area_m2) / self.tau_a
    /Users/oberrauch/oggm-fork/oggm/core/vascaling.py:1737: RuntimeWarning: invalid value encountered in double_scalars
      self.tau_a = self.tau_l * self.area_m2 / self.length_m ** 2
    /Users/oberrauch/oggm-fork/oggm/core/vascaling.py:115: RuntimeWarning: invalid value encountered in true_divide
      / (temp_grad * (max_hgt - min_hgt)))
    /Users/oberrauch/oggm-fork/oggm/core/vascaling.py:123: RuntimeWarning: invalid value encountered in double_scalars
      prcp_solid *= (1 + prcp_grad * (mean_hgt - ref_hgt)) * f_solid
    ```

  - RGI60-11.00757

    ```
    /Users/oberrauch/oggm-fork/oggm/core/vascaling.py:1781: RuntimeWarning: invalid value encountered in double_scalars
      self.dA = ((self.volume_m3 / self.ca) ** (1 / self.gamma)
    /Users/oberrauch/oggm-fork/oggm/core/vascaling.py:1786: RuntimeWarning: invalid value encountered in double_scalars
      self.dL = ((self.volume_m3 / self.cl) ** (1 / self.ql)
    ```

- Random climate with positive temperature bias

- - RGI60-11.00674
  - RGI60-11.00675
  - RGI60-11.00690
  - RGI60-11.00702
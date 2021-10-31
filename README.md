# MSc Data Science and Analytics (University of Leeds)
# Dissertation Module: MATH5872M


Topic: <b> Sample Size Re-estimation in Randomised trials </b>
<br>
<br>

This repository is a collection of all the codes, collected data and plots for the dissertation submitted in accordance with the master's degree. 
Please follow the below table for the structure of the files stored in this repository. 

<br>
<br>
<br>


<table>
<thead>
  <tr>
    <th>Folder</th>
    <th>Sub-folder</th>
    <th>File name</th>
    <th>Additional Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="3">Allocation Ratio Analysis</td>
    <td>-</td>
    <td>allocation_ratio_analysis.R</td>
    <td>R code for analyzing effects of change in allocation ratio on sample size re-estimation</td>
  </tr>
  <tr>
    <td>data</td>
    <td>multiple files</td>
    <td>Data files with data collected during simulations for different values of allocation ratio (k)</td>
  </tr>
  <tr>
    <td>plots</td>
    <td>multiple files</td>
    <td>Confidence interval plots using the collected data for allocation ratio analysis</td>
  </tr>
  <tr>
    <td rowspan="3">Method 1 - Chow et al. (2017)</td>
    <td>-</td>
    <td>chow_wang_shao.R</td>
    <td>R code for simulating the method proposed by the authors to analyze the method and collect relevant data</td>
  </tr>
  <tr>
    <td>data</td>
    <td>multiple files</td>
    <td>Data files with data collected during simulations based on the set parameters</td>
  </tr>
  <tr>
    <td>plots</td>
    <td>multiple files</td>
    <td>Confidence interval plots using the collected data for analysis of the method proposed by the authors</td>
  </tr>
  <tr>
    <td rowspan="4">Method 2 - Betensky and Tierney (1997)</td>
    <td>-</td>
    <td>betensky_tierney_M_value_analysis.R</td>
    <td>R code for simulating the method proposed by the authors to study the parameter 'M'</td>
  </tr>
  <tr>
    <td>-</td>
    <td>betensky_tierney.R</td>
    <td>R code for simulating the method proposed by the authors to analyze the method and collect relevant data</td>
  </tr>
  <tr>
    <td>data</td>
    <td>multiple files</td>
    <td>Data files with data collected during simulations based on the set parameters</td>
  </tr>
  <tr>
    <td>plots</td>
    <td>multiple files</td>
    <td>Confidence interval plots using the collected data for analysis of the method proposed by the authors</td>
  </tr>
  <tr>
    <td rowspan="3">Method 3 - Bristol and Shurzinske (2001)</td>
    <td>-</td>
    <td>bristol_shurzinske.R</td>
    <td>R code for simulating the method proposed by the authors to analyze the method and collect relevant data</td>
  </tr>
  <tr>
    <td>data</td>
    <td>multiple files</td>
    <td>Data files with data collected during simulations based on the set parameters</td>
  </tr>
  <tr>
    <td>plots</td>
    <td>multiple files</td>
    <td>Confidence interval plots using the collected data for analysis of the method proposed by the authors</td>
  </tr>
  <tr>
    <td rowspan="3">Method 4 - Kieser and Friede (2003)</td>
    <td>-</td>
    <td>bristol_shurzinske.R</td>
    <td>R code for simulating the method proposed by the authors to analyze the method and collect relevant data</td>
  </tr>
  <tr>
    <td>data</td>
    <td>multiple files</td>
    <td>Data files with data collected during simulations based on the set parameters</td>
  </tr>
  <tr>
    <td>plots</td>
    <td>multiple files</td>
    <td>Confidence interval plots using the collected data for analysis of the method proposed by the authors</td>
  </tr>
  <tr>
    <td>Visualization</td>
    <td>-</td>
    <td>visualizations.R</td>
    <td>R code for all the created plots in the analysis</td>
  </tr>
</tbody>
</table>






<br>
<br>
<br>
<br>
<br>
<br>


<b>References</b>


Chow S.,  Shao J.,  Wang H.  and Lokhnygina Y. (2017). Sample size calculations in clinicalresearch. <i>Chapman and Hall/CRC</i>.
  
Betensky A. and Tierney C. (1997). An examination of methods for sample size recalculation during an experiment. <i>Statistics in Medicine</i>, 16(22), 2587–2598.

Bristol D. and Shurzinske L. (2001). Blinded Sample Size Adjustment. <i>Drug Information Journal</i>. 35. 1123-1130.

Kieser M. and Friede T. (2003). Simple procedures for blinded sample size adjustment that do not affect the type I error rate. <i>Statistics in Medicine</i>. 15;22(23): 3571-81.


## Overview

A survey was conducted to rank schools given some data about the students. Our goal is to rank the schools based on the quality of teaching and learning that happens at each school, rather than ranking schools with an influence on school funding. For the analysis, we will use two different JAGS models to fit the "student" data and the "school" data. The two JAGS models are both mixed-effects models with either a student covariate, or a student and a school covariate.

## Variable Descriptions

- `Y:` The response.
- `LRT:` LRT.
- `VR1:` VR1.
- `VR2:` VR2.
- `Gender:` The gender of the individual.
- `School:` The school that the individual attends.
- `CE:` CE.
- `RC:` RC.
- `Other:` Other.
- `Girls:` Whether the individual is a girl.
- `Boys:` Whether the individual is a boy.

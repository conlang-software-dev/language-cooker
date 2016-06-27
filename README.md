# language-cooker

Automatic generation of typologically-plausible language descriptions using R.

Based on the morphosyntax features coded in the [World Phonotactics Database](http://phonotactics.anu.edu.au/qtp/?db=wals-latest) (thanks to [Mark Donohue](https://anu-au.academia.edu/MarkDonohue) for making a download of the data available); see [WALS](http://wals.info/chapter) for descriptions of features and feature values.

The program works by computing a correlation matrix among all the features, generating a vector of uniformly-distributed random numbers that follow the same correlations, and comparing each element of the vector with the prevalence of the corresponding feature in the relevant sample.

It is possible to generate a language that is typologically representative (i.e. based on the feature distributions in a typologically balanced sample); globally representative (based on global feature distributions without sampling); or representative of a particular language family, macro-area (~ continent), or country.

More concretely, you can generate languages that look like they might have come from anywhere on the planet; or you can generate a language that “looks Uralic”, or “looks African”, or “looks Nepalese”.

There is currently a demo of the program up on [DataJoy](https://www.getdatajoy.com/project/5766150870b7197577383711); see the README file there to get started. (Note: The program will from time to time produce feature combinations that seem logically impossible. Use common sense. Also, don't pay too much attention to the “Your language resembles X” statements; the generated language will usually be most similar to the few most average languages in the sample, but this is because most languages are _quite far_ from the average.)

The name “Language Cooker” is inspired by [Type Cooker](http://www.typecooker.com), a website that randomly generates descriptions of typefaces or typeface families, to provide exercises for students of typeface design. Likewise, the Language Cooker could be used to produce exercises for students of language design, i.e. conlangers.

See also “[Linguistx](http://linguistx.tumblr.com)”’s [Colorless](http://linguistx.tumblr.com/colorless), a language generator that generates phonologies as well (based on correlations in [UPSID](http://web.phonetik.uni-frankfurt.de/upsid.html)), and converts the language descriptions into natural language rather than expressing them as lists of numerical features. (However, there is no explanation of how the generator works, or how the correlations were arrived at.) See also Alex Fink’s [Gleb](https://github.com/alexfink/random_language), which takes a sophisticated approach to procedurally generating phonologies (which I don’t understand yet).

# **Greetings**

Good morning from me too, and thanks for joining today for the presentation and defense of my Masters Thesis titled "Testing the importance of explicit glacier dynamics for mountain glacier change projections".

# **Table of contents**

I'll start with a brief introduction,  providing you an outlook on what to expect from this presentation, talking about the motivation and goals of my thesis, before going into some more details about the glacier models used hereafter.

Then I'll present you the main results of this work, starting with a single glacier test case, continuing with regional application, and a sensitivity analysis and ending with 21st century projection of glacier volume change for the Alps and High Mountain Asia.

I'll conclude the presentation with answering the research questions presented in the introduction.

# **Introduction**

Let me start with the the most obvious question: WHAT did I do?!

This work is  basically a miniature model intercomparison study, comparing the **volume/area scaling model**

originally presented by **Marzeion** and colleagues in 2012 used to estimate past and future sea level rise contribution of glaciers

to the **OPEN GLOBAL GLACIER MODEL**.

The OGGM is a flowline model and was developed out of the same paper/same idea but has  a dedicated ice physics module. But, more on that later.

How did I do it?!

Re-implementing the VAS model in the OGGM framework **allows for controlled boundary conditions**. That means that both models use the **same input data** like climate files, glacier outline, digital elevation models, the same mass balance model and calibration strategy. Thereby it is assured that **all differences in model output** will arise from differences in the dynamic response.

And Why did I do it?!

I'm sure everybody here today is aware of the **past and ongoing climate change and its implications**. Amongst others, the rise in **average global temperature** results in an increased ice melt of glaciers and ice caps, which in turn contribute to the global sea level rise.

While sea level rise may be the overarching motivation, in this work we also wanted to **asses strength and weaknesses of both models**.

The increase in computational power and data availability over the last years has allowed for new and advanced models like the OGGM. I'm investigate what **additional information is added by increasing** the model complexity.

In summary I tried to **answer the following two questions**:

- *What are the differences in dynamic response between the volume/area scaling model and the flowline model?*

- *How do those differences influence ice volume projections for the 21st century for mountain regions?*

So, without further ado, let's step into the definition of the volume area scaling model:



# **Model structure**

There are two major steps in the workflow. Before the actual forward run, the **model glacier model must be initialized**. So we'll look at this initialization process first:

Each glacier is **modeled individually**, here I'll use the Hintereisferner as example.

The model starts with the glacier **outlines and a digital elevation model**  from which it calculates the initial glacier surface area.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210413231135803.png" alt="image-20210413231135803" style="zoom:25%;" />

At that point, the **scaling relations come into play**. The volume/area scaling is a **simple power law relation** with a scaling constant and a scaling exponent.

If we plug in the initial area we get an estimate for the initial ice volume.

There is an **analog scaling relation between glacier length and ice volume**, with different scaling parameters.

If we **invert** this relation and plug in the above computed initial ice volume, we can **solve for the initial glacier length.**

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210413231309949.png" alt="image-20210413231309949" style="zoom:25%;" />

And so that's how the **model represents a glacier**, with its initial values for area, volume and length.

Additionally, we have information about the **maximum surface elevation**, which is assumed to be constant over the lifetime of the glacier, and the elevation of the **glacier tongue, the so called terminus elevation**, which can change over the course of time and is here denoted with 0 for its initial value.

The **elevation range in combination with** the glacier lengths informs about the **average surface slope**.

# **Forward run**

Now let's step into the forward run, which basically consists of the **interplay between mass balance model and ice dynamics model**.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210413231641230.png" alt="image-20210413231641230" style="zoom:25%;" />

Even though it is a **bit more complicated than that**.

We have to **update the glacier geometries** from one time step to another, so from this year t to the next year t+1.

We do so by **starting with computing the volume change** and the new volume, for which we **need the specific mass balance** first.

A glaciers mass balance is defined as the **difference between accumulation**, or mass gain, and **ablation**, or mass loss. The implemented mass balance model considers only accumulation in form of solid precipitation, so snow fall, and ablation only due to ice melt.

The **mass balance equation** for the specific mass balance computes the **difference between solid precipitation and melted ice** , whereby the amount of ice melt is a **function of positive melting temperature and a glacier specific temperature sensitivity.** This difference is summed over **all month of a hydrological year** and corrected by a potential mass balance bias.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210413231919265.png" alt="image-20210413231919265" style="zoom:25%;" />

These two **glacier specific parameters** must be calibrated. The calibration is specific to this mass balance model, but since both the scaling model and the flowline model use it, I wont go into **detail now but again am happy to do so after the presentation**.

The specific mass balance is a **glacier wide average value per unit area**, hence it must be multiplied with the surface area and scaled with the ice density to get a value for **total volume change**.

What is important is that the **volume change happens instantenously**, so from one time step to the next. This means that next years volume is just **the sum of this years volume and the computed volume change**, which can be positive or negative depending on the specific mass balance

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210413232213910.png" alt="image-20210413232213910" style="zoom:25%;" />

However, **changes in overall volume drive changes in glacier length and surface area.**

So from the **newly computed volume** we can derive changes in length and surface area, but for that we **need the response time scales.**

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210413232451928.png" alt="image-20210413232451928" style="zoom:25%;" />

If we would just plug the new volume in the volume/area scaling relation we would get the **equilibrium area**, so the new area the **glacier reaches eventually** after it fully adjusted to the new volume. However, it **takes time** for those changes to take place.

Hence, the **area changes follow response time scaling**, another scaling relation that **must be used in conjuncture** with volume/area scaling **when applied to transient glaciers**.

So this **actual area change** is the **total possible area change**, so the difference between the equilibrium area and the **current area**, **scaled by the area response time**. The new area is then again just the **sum of current area and computed area change.**

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210413232917382.png" alt="image-20210413232917382" style="zoom:25%;" />

The same holds true for **glacier length**, with the corresponding scaling parameters and **the length response time.**

**Last but not least**, changes in glacier length then **inform changes in terminus elevation**.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210413233032746.png" alt="image-20210413233032746" style="zoom:25%;" />

The terminus elevation is a **linear function of the glacier length**, so assuming an average gradient of 10% and a length change of 100m the terminus elevation would decrease by 10 meters. 

# **Flowline model**

Now **briefly** to the **Open Global Glacier Model**.

The OGGM is a **flowline model**, it represents the glacier with a main flowline and **possible tributaries**.

It is **1.5 dimensional**, that means that while the **ice flows only along one direction** namely along the flowline, each grid point on the flowline has **information about the glacier width, ice thickness and bed shape**. This makes the model **geometry aware.**

It uses a more complex ice thickness inversion scheme, which also allows to compute **distributed ice thickness estimation.**

The ice dynamics module **is also  more advanced**, solving the shallow ice equations on a staggered grid for the flowline glacier.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414083219812.png" alt="image-20210414083219812" style="zoom:25%;" />

Here we have the **flowlines with the transectional area** of one grid point.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414083241103.png" alt="image-20210414083241103" style="zoom:25%;" />

The volume conservation equation equates the **change in cross sectional area** of this transect, to the **difference of specific point mass balance** and **ice flux divergence**. Due to time constraints I'll spare you the rest of the equations, but again I'm happy to elaborate on them afterwards.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414083319515.png" alt="image-20210414083319515" style="zoom:25%;" />

The equations are solved by the **FluxBasedModel. As the name might suggest**, it solves for the flux divergence term in the equation, which **allows to use the same equation** and solver for **all three different types of bed shapes.**

This is a **very superficial description** of the OGGM flowline model, but I hope you got an idea about the **added complexities in comparison to the scaling model.**

# **Single glacier test case**

Now let's step into the results:

The first experiment I'd like to show you is **an initial test case**, here again on the example of Hintereisferner.

While scaling relations should be **applied to collections of glaciers** and not single glaciers, this single glacier test case **nicely highlights certain characteristics** of the scaling model.

We ran **equilibrium experiments** with three different temperature biases. What does that mean?!

We start with a glacier in equilibrium and the **corresponding climate, which is then perturbed** with a step change in temperature. Here for the **cooling scenario we decrease** the temperature by 0.5°C while keeping everything else the same. This **in turn results in a glacier growth** until a new equilibrium is reached, eventually.

The same is done with a positive temperature bias of 0.5°C, which results in a glacier retreat.

**Additionally, we have an equilibrium scenario as a sort of control run.**

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414083541466.png" alt="image-20210414083541466" style="zoom:25%;" />

This is done using two different mass balance models: the constant mass balance model and the random mass balance model. As the names might suggest, the former keeps the climate constant as depicted here, while the later introduces a natural random year to year variability. This is obviously more realistic, but also results in a noisier output.

# **Hintereisferner - ice volume**

This figure has a bit much information, so let me unpack it for you.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414083733275.png" alt="image-20210414083733275" style="zoom:25%;" />

Here we **have the temporal evolution of ice volume for Hintereisferner**, modeled by the **volume/area scaling** model on the left and by the **flowline model** on the right. The **x-axis shows times**, while the y-axis shows the **relative ice volume, normalized** with the initial value at the **start of the simulation** for better comparability between both models.

Then we have the results of **constant climate scenarios in dashed lines** and the results of the **random climate scenarios in solid lines**.

If we look at the random climate scenarios we see that there are some **differences in the models response to natural climate variability**, more low frequency variability for the flowline mode and more high frequency variability for the scaling model.

We will take a **closer look at that later**, but let's ignore it **for now and focus on the constant** climate scenarios.

There are three things I want to draw your attention to:

1. **Volume changes** estimated by the scaling model are much smaller than the ones estimated by the flowline model.

   The scaling model estimates a change of +/- 20%, while the flowline model estimates an increase in ice volume of over 70% and a loss of ice of over 40%.

2. This **brings as to the second point**. The scaling model produces **highly symmetrical results, due the missing mass balance** elevation feedback. If this feedback is turned off for the flowline model, the results looks similarly symmetric.

3. And finally, while the flowline model shows **an exponential or rather asymptotical changes in ice volume**, the scaling model **shows an oscillating behavior** comparable to a damped harmonic oscillator, which is hardly a physical result.



The same basic observations hold true for the **surface area, though the oscillation is more pronounced**

and the **glacier length, though the difference in estimated changes is even more drastic** even more so for the **absolute values**. The glacier length of the scaling model should be considered more as a model parameters rather than an actual physical value.

**Nevertheless**, now we are going to take a closer look at the **glacier models response to natural climate variability** using the **length changes** under random climate scenarios as input.



# **Hintereisferner - natural length fluctuation**

This figure shows **length anomalies with respect to the corresponding equilibrium value** for a period of 4000 years of out of the total of 23 000 years of simulation.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414084411476.png" alt="image-20210414084411476" style="zoom:25%;" />

**Each panel represents one temperature bias**, same as before we have a **cooling** scenario, a **equilibrium** scenario and a **warming** scenario. However, the initial **3 000 year in which the glaciers adjusted** to the new climate were chopped, hence the **temperature biases can be seen as labels for glacier size**. So we have a big, normal and small size of the same glacier.

Then again, we have the results of the **volume area scaling** model in yellow, orang and red lines and on the left y-axis.

And the results of the **flowline model in blueish lines and on the right y-axis**.

**Note the differences in scale**, the length **anomalies estimated by the flowline model are more than twice** that of the scaling model.

**Two more things**

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414085107625.png" alt="image-20210414085107625" style="zoom:25%;" />

To quantify this we perform a **power spectral density analysis of the glacier length signal**. The power spectral density is related to the **discrete Fourier transformation** and used characterize random processes. It describes how the **power of a signal is distributed over frequency**.

We see **higher and rather constant power for frequencies** up to to about here (3 · 10−3), which corresponds to signals with a **period of over 300 years**. After that, the power **decreases with increasing frequency**.

The **power spectral density of the corresponding mass balance**, however,  we see a **constant line in the log/log space**. This means that the power is **uniformly distributed, no dominant frequencies, no periodicity**. So we have a random input signal, and a non uniform output signal.

This makes **intuitive sense**, since changes in glacier length are **mainly driven by long term climatic trends and less by inter-annual variabilities** in the climatic forcing. This so **called low-pass filter behavior is characteristic** for glaciers and both models are able to represent it. A low pass filter passes low frequencies, while it block or attenuates higher frequencies.

However, we once again **see no difference in the different runs of the scaling model**. While the flowline model allows the differently sized glaciers to evolve differently.

# **Alpine runs**

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414085255187.png" alt="image-20210414085255187" style="zoom:25%;" />

As I said before, scaling relations should be **applied to collections of glaciers at least on a regional scale** to provide more than just **order of magnitude estimations**. Hence, I repeated the **equilibrium experiments for the entirety of all Alpine glaciers**.

Here we see the **aggregate ice volume** of the Alps, same as before for Hintereisferner. All observations made for the single glacier still hold true:

- the **scaling model estimated less change** in ice volume than the flowline model
- it produces **symmetrical** results
- and it shows the **oscillating** behavior.

We see that on a regional scale there is little difference between the random and constant climate scenario.

**To its defense**, the **initial ice volume estimation of the scaling model fits quite good** with the latest consensus ice volume estimation published Farinotti and colleagues in 2019. The **flowline model in its out of the box state tends to overestimate** the ice volume, but at least it can be **calibrated**.

Hence, in the next slide we look **at a sensitivity analysis and possible tuning parameters** for the scaling model.

For that I used only the **warming scenario of the alpine run and varied the key parameters of both scaling relations.**

# **Sensitivity experiments**

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414085930396.png" alt="image-20210414085930396" style="zoom:25%;" />

Let's first look at the **sensitivity to the response time scaling**.

The plot shows the **relative aggregate ice volume for all alpine glaciers under a warming scenario of 0.5 degree celsius**. Thereby I scaled the **model internal time scale with a simple linear factor**.

It is apparent, that the **time scales act as damping factor for the oscillation**. **Halving** the time scale, which means letting the **glacier react faster**, weakens the **oscillation and vice versa**.

The **equilibrium values are not affected**, for that we have to look at the volume/area and volume length scaling.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414090149133.png" alt="image-20210414090149133" style="zoom:25%;" />

The plot **again shows three different runs:**

- one with the same **global constants** as used for all runs up until now
- one with **custom scaling constants** while keeping the scaling exponent fixed
- and one with **custom constants and exponents**.

Changing the constants only results in **different absolute values**, which in **turn slightly affects the time scales** and thereby the temporal evolution. However, the **relative values normalized with the initial volume** do not changes.

To do so, we have to use **custom scaling constants and exponents**. However, the **differences in equilibrium volume are relatively minor** and still **far off from the one estimated by the flowline** model.

This can be **seen both way:**

- while there is **some potential for model calibration** it is unlikely to get **comparable (absolute) results** for the scaling model and the flowline model.
- however, the model is **quite insensitive to the scaling parameters**.

# **Projections**

Finally, after all the **idealized experiments** we arrive at the actual projections for the current century.

This was done by **using CRU time series climate data** for the run up to 2020, after which we used climate data from **fifteen general circulation model** of the newest **climate model intercomparison project**, with four different **SHARED SOCIO ECONOMIC PATHWAYS** each. which is just a fancy word for saying scenarios.

To **further increase the accuracy of the results** and assure similar region wide surface mass balance estimates between both models, the **mass balance residual of both models is corrected to match the regional geodetic observations** presented by Hugonnet et al. (2020).

The initial  idea was to limit this study to the Alps, but then we decided to include other mountainous regions with more and bigger glaciers.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414090505031.png" alt="image-20210414090505031" style="zoom:25%;" />

First we look at the results for Central Europe, which besides the Alps contains some glaciers in the Apennine and the Pyrenees, again the scaling model on the left and the flowline model on the right.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414090550573.png" alt="image-20210414090550573" style="zoom:25%;" />

The volume is **shown as fraction of the 2020 volume** for both models, but **the secondary y-axis shows the corresponding absolute values in units of sea level rise milli meters**. Those axes **differ between both models since the initial values** in 2020 are different.

The **different SSPs are denoted in different colors**, whereby the **shading represents the uncertainty bounds**, so about 68% of all estimations lie inside the shaded area.

I'd like to point out that despite  **basically identical mass balance estimations**, the flowline model **glaciers loose more ice mass** (relative and absolute) over the first 20 years than the volume/area scaling model glaciers

After **2020 the ice loss continues with the same constant rate** before attenuating in the second half of the 21st century. For **all but the lowest forcing scenarios, all or almost all ice will be gone** by the end of this century.

The ensemble averages equal to a sea level rise between **0.16 mm SLE** for SSP1-2.6 and and **0.24mm SLE SSP5-8.5** 

The **scaling model predicts much less relative ice loss**, especially for the SSPs with lower forcings.

so much so that the **final estimates of one model are outside the margin of error of the other**

**Only for SSP5-8.5** does the estimated absolute mass loss of **0.21mmSLE** compare to the flowline model, while for **SSP1-2.6 with 0.08mmSLE** it is only half of the flowline model estimate.

<img src="/Users/oberrauch/Library/Application Support/typora-user-images/image-20210414091758293.png" alt="image-20210414091758293" style="zoom:25%;" />

If we look at Central Asia, which has **almost 15 times as many glaciers** as central Europe and especially much bigger ones, we see a different picture.

This results in **less loss in relative terms but much more in absolute units**. Also the overall evolution is different. At the **end of the century, there is still about half** of the ice left and hence the glaciers are on a downward trajectory.

Again, the flowline model predicts more ice loss (before and after 2020) than the volume/area scaling model.

While the **relative values can be drastically different**, the **absolute values in SLE are in better agreement,** to the end that **uncertainty estimates of both models overlap** for most regions and SSPs 
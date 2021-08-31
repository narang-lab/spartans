import React from 'react';
import clsx from 'clsx';
import styles from './HomepageFeatures.module.css';

const FeatureList = [
  {
    title: 'Spatially-Resolved, Momentum-Exact',
    Svg: require('../../static/img/feature_01.svg').default,
    description: (
      <>
        SpaRTaNS solves the steady-state Boltzmann transport equation in 
        three real-space and three momentum-space dimensions.
      </>
    ),
  },
  {
    title: 'General Framework',
    Svg: require('../../static/img/feature_02.svg').default,
    description: (
      <>
        SpaRTaNS is carrier-agnostic with applications in fields like
        phonon transport and <q>hydrodynamic</q> electron transport.
      </>
    ),
  },
  {
    title: 'Open Source & Parallel',
    Svg: require('../../static/img/feature_03.svg').default,
    description: (
      <>
        SpaRTanS makes use of multiple processors using the MPI framework,
        and is embarassingly-parallel across real-space tetrahedra.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Svg className={styles.featureSvg} alt={title} />
      </div>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}

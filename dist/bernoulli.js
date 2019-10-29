/*
  Author:  David Faunce
  Date:    10/27/2019
  Purpose: Solves Bernoulli's Fluid Function for incompressible fluid flow
  
  Input:   Inlet Properties (Pressure, Velocity, Height, Fluid Density), 
           Outlet Properties (same properties as inlet)
  Output:  The solver finds the missing data for the Inlet and Outlet 
           properties and returns an object of all the properties.
*/
$(function() {

  function bernoulli(options) {
    var $o = $.extend({
      returnSIUnits: true,
      sameMedium: true,
      neglectHeight: true,
      g: {val: 9.81, unit:"m/s2"},
      p1: {val: null, unit:null},
      rho1: {val: null, unit:null},
      v1:{val:null, unit:null},
      Q1: {val: null, unit: null},
      A1: {val: null, unit: null},
      h1: {val:null, unit:null},
      p2: {val: null, unit:null},
      rho2: {val: null, unit:null},
      v2: {val:null, unit:null},
      Q2: {val: null, unit: null},
      A2: {val: null, unit: null},
      h2: {val:null, unit:null},
      Kin_1: null,
      Pot_1: null,
      Kin_2: null,
      Pot_2: null
    }, options);
    
    //convert to SI units
    $o.p1.val = ($o.p1.val != null) ? jsConvert($o.p1.val, $o.p1.unit, "P") : $o.p1.val;
    $o.p2.val = ($o.p2.val != null) ? jsConvert($o.p2.val, $o.p2.unit, "P") : $o.p2.val;
    
    $o.rho1.val = ($o.rho1.val != null) ? jsConvert($o.rho1.val, $o.rho1.unit, "kg/m3") : $o.rho1.val;
    $o.rho2.val = ($o.rho2.val != null) ? jsConvert($o.rho2.val, $o.rho2.unit, "kg/m3") : $o.rho2.val;
    
    $o.v1.val = ($o.v1.val != null) ? jsConvert($o.v1.val, $o.v1.unit, "m/s") : $o.v1.val;
    $o.v2.val = ($o.v2.val != null) ? jsConvert($o.v2.val, $o.v2.unit, "m/s") : $o.v2.val;
    
    $o.Q1.val = ($o.Q1.val != null) ? jsConvert($o.Q1.val, $o.Q1.unit, "m3/s") : $o.Q1.val;
    $o.Q2.val = ($o.Q2.val != null) ? jsConvert($o.Q2.val, $o.Q2.unit, "m3/s") : $o.Q2.val;
    
    $o.A1.val = ($o.A1.val != null) ? jsConvert($o.A1.val, $o.A1.unit, "m2") : $o.A1.val;
    $o.A2.val = ($o.A2.val != null) ? jsConvert($o.A2.val, $o.A2.unit, "m2") : $o.A2.val;
    
    $o.h1.val = ($o.h1.val != null) ? jsConvert($o.h1.val, $o.h1.unit, "m2") : $o.h1.val;
    $o.h2.val = ($o.h2.val != null) ? jsConvert($o.h2.val, $o.h2.unit, "m2") : $o.h2.val;
    
    if ($o.neglectHeight) {
      if ($o.h1.val == null && $o.h2.val != null) {
        $o.h1.val = $o.h2.val;
      }
      else if ($o.h1.val != null && $o.h2.val == null) {
        $o.h2.val = $o.h1.val;
      }
      else if ($o.h1.val == null && $o.h2.val == null) {
        $o.h1.val = 0;
        $o.h2.val = 0;       
      }
      
      $o.h1.unit = "m";
      $o.h2.unit = "m";
    }
    
           
    function setRho() {
      if ($o.sameMedium && ($o.rho1.val == null || $o.rho2.val == null)) {
       $o.rho1.val = ($o.rho2.val || $o.rho1.val);
       $o.rho2.val = ($o.rho1.val || $o.rho2.val);
      } 
    }
    
    setRho();
   
    //Solve Continuity
    function solveContinuity(cnt) {
      if (cnt >= 3) {
        return;
      }      
      if ($o.Q1.val == null) {
        if ($o.Q2.val != null) {
          $o.Q1.val = $o.Q2.val;
        }
        else if ($o.A1 != null && $o.v1 != null) {
          $o.Q1.val = $o.A1.val * $o.v1.val;
        }        
      }
      
      if ($o.Q2.val == null) {        
        if ($o.Q1.val != null) {
          $o.Q2.val = $o.Q1.val;
        }
        else if ($o.A2 != null && $o.v2 != null) {
          $o.Q2.val =  $o.A2.val * $o.v2.val;
        }
      }   
      
      $o.A1.val = ($o.A1.val == null && $o.Q1.val != null && $o.v1.val != null) ? $o.Q1.val / $o.v1.val : $o.A1.val;
      $o.v1.val = ($o.v1.val == null && $o.Q1.val != null && $o.A1.val != null) ? $o.Q1.val / $o.A1.val : $o.v1.val;
      
      $o.A2.val = ($o.A2.val == null && $o.Q2.val != null && $o.v2.val != null) ? $o.Q2.val / $o.v2.val : $o.A2.val;
      $o.v2.val = ($o.v2.val == null && $o.Q2.val != null && $o.A2.val != null) ? $o.Q2.val / $o.A2.val : $o.v2.val;
      
      cnt++;
      solveContinuity(cnt);
    }
    
    solveContinuity(0);
    
    //Solve Kinetic Energy of Location 1
    function solveKin_1() {
      if ($o.Kin_1 == null && $o.rho1.val != null && $o.v1.val != null) {
        $o.Kin_1 = 0.5 * $o.rho1.val * $o.v1.val * $o.v1.val;
      }
      return $o.Kin_1 === true;
    }
    
    //Solve Kinetic Energy of Location 2
    function solveKin_2() {
      if ($o.Kin_2 == null && $o.rho2.val != null && $o.v2.val != null) {
        $o.Kin_2 = 0.5 * $o.rho2.val * $o.v2.val * $o.v2.val;
      }
      return $o.Kin_2 === true;
    }
    
    //Solve Potential Energy of Location 1
    function solvePot_1() {
      if ($o.Pot_1 == null && $o.rho1.val != null && $o.g.val != null && $o.h1.val != null) {
        $o.Pot_1 = $o.rho1.val * $o.g.val & o.h1.val;
      }
      return $o.Pot_1 === true;
    }
    
    //Solve Potential Energy of Location 2
    function solvePot_2() {
      if ($o.Pot_2 == null && $o.rho2.val != null && $o.g.val != null && $o.h2.val != null) {
        $o.Pot_2 = $o.rho2.val * $o.g.val & o.h2.val;
      }
      return $o.Pot_2 === true;
    }
    
    function solveEnergies() {
      return (solveKin_1 && solveKin_2 && solvePot_1 && solvePot_2);
    }
    
    function solveP1() {
      if ($o.p1.val == null) {
        if (solveEnergies() && $o.p2.val != null) {
          $o.p1.val = $o.p2.val + $o.Kin_2 + $o.Pot_2 - $o.Kin_1 - $o.Pot_1;
        }        
      }
    }
    function solveP2() {
      if ($o.p2.val == null) {
        if (solveEnergies() && $o.p1.val != null) {
          $o.p2.val = $o.p1.val + $o.Kin_1 + $o.Pot_1 - $o.Kin_2 - $o.Pot_2; 
        }
      }
    }
    function solveH1() {
      if ($o.h1.val == null) {
        solveEnergies();
        if ($o.p1.val != null && $o.p2.val != null && $o.Kin_1 != null && $o.Kin_2 != null && $o.Pot_2 != null && $o.rho1.val != null) {
          $o.h1.val = (($o.p2.val - $o.p1.val) + $o.Kin_2 - $o.Kin_1 + $o.Pot_2) / ($o.rho1.val * $o.g.val);
        }
      }
    }
    function solveH2() {
      if ($o.h2.val == null) {
        solveEnergies();
        if ($o.p1.val != null && $o.p2.val != null && $o.Kin_1 != null && $o.Kin_2 != null && $o.Pot_1 != null && $o.rho2.val != null) {
          $o.h1.val = (($o.p1.val - $o.p2.val) + $o.Kin_1 - $o.Kin_2 + $o.Pot_1) / ($o.rho2.val * $o.g.val);
        }
      }
    }
    
    function solveV1() {
      if ($o.v1.val == null) {
        solveEnergies();
        if ($o.p1.val != null && $o.p2.val != null && $o.Kin_2 != null && $o.Pot_1 != null && $o.Pot_2 != null && $o.rho1.val != null) {
          var s = 2 * (($o.p2.val - $o.p1.val) + $o.Kin_2 + $o.Pot_2 - $o.Pot_1) / $o.rho1.val;
          $o.v1.val = Math.sqrt(s);
        }
      }
    }
    
    function solveV2() {
      if ($o.v2.val == null) {
        solveEnergies();
        if ($o.p1.val != null && $o.p2.val != null && $o.Kin_1 != null && $o.Pot_1 != null && $o.Pot_2 != null && $o.rho2.val != null) {
          var s = 2 * (($o.p1.val - $o.p2.val) + $o.Kin_1 + $o.Pot_1 - $o.Pot_2) / $o.rho2.val;
          $o.v2.val = Math.sqrt(s);
        }
      }
    }
    
    function solveRho2() {
      if ($o.rho2.val == null) {
        solveEnergies();
        if ($o.p1.val != null && $o.p2.val != null && $o.Kin_1 != null && $o.Pot_1 != null && $o.v2.val != null && $o.h2.val != null) {
          $o.rho2.val = (2 * ($o.p1.val - $o.p2.val) + $o.Kin_1 + $o.Pot_1) / (($o.v2.val * $o.v2.val) + (2 * $o.g.val * $o.h2.val));
        }
      }
    }
    function solveRho1() {
      if ($o.rho1.val == null) {
        solveEnergies();
        if ($o.p1.val != null && $o.p2.val != null && $o.Kin_2 != null && $o.Pot_2 != null && $o.v1.val != null && $o.h1.val != null) {
          $o.rho1.val = (2 * ($o.p2.val - $o.p1.val) + $o.Kin_2 + $o.Pot_2) / (($o.v1.val * $o.v1.val) + (2 * $o.g.val * $o.h1.val));
        }
      }
    }
    
    
    function init(cnt) {      
      if (cnt >= 5) {
        return $o;
      }
      solveContinuity(0);
      solveEnergies();
      solveP1();
      solveP2();
      solveV1();
      solveV2();
      solveH1();
      solveH2();
      solveRho1();
      solveRho2();
      cnt++;
      init(cnt);
    }   
    
    init(0);
  }
  
  
  function foo() {
    var obj = {
      returnSIunits: true,
      sameMedium: true,
      neglectHeight: true,
      g: {val: 9.81, unit:"m/s2"},
    
      p1: {val: null, unit:"P"},
      rho1: {val: 1000, unit:"kg/m3"},
      v1:{val:1.96, unit:"m/s"},
      Q1: {val: null, unit: null},
      A1: {val: null, unit: null},
      h1: {val:0, unit:"m"},
  
  
      p2: {val: 101000, unit:"P"},
      rho2: {val: null, unit:"kg/m3"},
      v2: {val:25.5, unit:"m/s"},
      Q2: {val: null, unit: null},
      A2: {val: null, unit: null},
      h2: {val:0, unit:"m"} 
    };
    console.log(bernoulli(obj));
  }
  
  foo();
  
});

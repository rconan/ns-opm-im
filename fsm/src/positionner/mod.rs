use dosio::{io::Tags, ios, DOSIOSError, Dos, IOTags, IOVec, IO};
use simulink_binder::import;

import! {M2_POS_Control,
r##"	 
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
/*
 * File: M2_POS_Control.h
 *
 * Code generated for Simulink model 'M2_POS_Control'.
 *
 * Model version                  : 1.900
 * Simulink Coder version         : 9.0 (R2018b) 24-May-2018
 * C/C++ source code generated on : Tue Aug 24 14:23:00 2021
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_M2_POS_Control_h_
#define RTW_HEADER_M2_POS_Control_h_
#include <string.h>
#include <stddef.h>
#ifndef M2_POS_Control_COMMON_INCLUDES_
# define M2_POS_Control_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                                 /* M2_POS_Control_COMMON_INCLUDES_ */

#include "M2_POS_Control_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T FBcontroller_states[126];     /* '<S1>/FB controller' */
} DW_M2_POS_Control_T;

/* Constant parameters (default storage) */
typedef struct {
  /* Expression: m2pos.Km2p_dec
   * Referenced by: '<S1>/Gain'
   */
  real_T Gain_Gain[1764];

  /* Expression: kron(eye(42),[1;-1])
   * Referenced by: '<S1>/kron(eye(42),[1;-1])'
   */
  real_T kroneye4211_Gain[3528];
} ConstP_M2_POS_Control_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T M2_pos_cmd[42];               /* '<Root>/M2_pos_cmd' */
  real_T M2_pos_FB[84];                /* '<Root>/M2_pos_FB' */
} ExtU_M2_POS_Control_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T M2_pos_act_F[84];             /* '<Root>/M2_pos_act_F' */
} ExtY_M2_POS_Control_T;

/* Real-time Model Data Structure */
struct tag_RTM_M2_POS_Control_T {
  const char_T * volatile errorStatus;
};

/* Block states (default storage) */
extern DW_M2_POS_Control_T M2_POS_Control_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_M2_POS_Control_T M2_POS_Control_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_M2_POS_Control_T M2_POS_Control_Y;

/* Constant parameters (default storage) */
extern const ConstP_M2_POS_Control_T M2_POS_Control_ConstP;

/* Model entry point functions */
extern void M2_POS_Control_initialize(void);
extern void M2_POS_Control_step(void);
extern void M2_POS_Control_terminate(void);

/* Real-time Model object */
extern RT_MODEL_M2_POS_Control_T *const M2_POS_Control_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S1>/m2_pos_en' : Eliminated nontunable gain of 1
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Note that this particular code originates from a subsystem build,
 * and has its own system numbers different from the parent model.
 * Refer to the system hierarchy for this subsystem below, and use the
 * MATLAB hilite_system command to trace the generated code back
 * to the parent model.  For example,
 *
 * hilite_system('ims_Build5pt1e/M2_POS_Control')    - opens subsystem ims_Build5pt1e/M2_POS_Control
 * hilite_system('ims_Build5pt1e/M2_POS_Control/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'ims_Build5pt1e'
 * '<S1>'   : 'ims_Build5pt1e/M2_POS_Control'
 */
#endif                                 /* RTW_HEADER_M2_POS_Control_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
"##}

impl<'a> IOTags for Controller<'a> {
    fn outputs_tags(&self) -> Vec<Tags> {
        vec![ios!(MCM2SmHexF)]
    }
    fn inputs_tags(&self) -> Vec<Tags> {
        ios!(M2poscmd, MCM2SmHexD)
    }
}
impl<'a> Dos for Controller<'a> {
    type Input = Vec<f64>;
    type Output = Vec<f64>;
    fn inputs(&mut self, data: Option<Vec<IO<Self::Input>>>) -> Result<&mut Self, DOSIOSError> {
        match data {
            Some(mut data) => {
                if let Some(
                    [IO::M2poscmd {
                        data: Some(m2_pos_cmd),
                    }, IO::MCM2SmHexD {
                        data: Some(m2_pos_fb),
                    }],
                ) =
                    <Vec<IO<Vec<f64>>> as IOVec>::pop_these(&mut data, ios!(M2poscmd, MCM2SmHexD))
                        .as_ref()
                        .map(|x| x.as_slice())
                {
                    for (k, &v) in m2_pos_cmd.into_iter().enumerate() {
                        self.m2_pos_cmd[k] = v;
                    }
                    for (k, &v) in m2_pos_fb.into_iter().enumerate() {
                        self.m2_pos_fb[k] = v;
                    }
                    Ok(self)
                } else {
                    Err(DOSIOSError::Inputs(
                        "FSM positioner M2poscmd and MCM2SmHexD not found".into(),
                    ))
                }
            }
            None => Err(DOSIOSError::Inputs(
                "None data passed to FSM positionner controller".into(),
            )),
        }
    }
    fn outputs(&mut self) -> Option<Vec<IO<Self::Output>>> {
        Some(vec![ios!(MCM2SmHexF(Vec::<f64>::from(&self.m2_pos_act_f)))])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positionner_feedback() {
        let mut ctrlr = Controller::new();
        let u = ios!(
            M2poscmd(vec![0f64; 42]),
            MCM2SmHexD({
                let mut v = vec![0f64; 84];
                v[0] = 1e-6;
                v[1] = -1e-6;
                v
            })
        );
        let y = ctrlr.in_step_out(Some(u));
        println!("y:\n{:#?}", y);
    }
}

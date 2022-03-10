stk.v.12.0
WrittenBy    STK_v12.0.0

BEGIN Chain

    Name		 LOS
    BEGIN Definition

        Object		 Place/Moscow
        Object		 Satellite/Satellite2
        Type		 Chain
        FromOperator		 Or
        FromOrder		 1
        ToOperator		 Or
        ToOrder		 1
        Recompute		 Yes
        IntervalType		 0
        ComputeIntervalStart		 0
        ComputeIntervalStop		 86400
        ComputeIntervalPtr		
        BEGIN EVENTINTERVAL
            BEGIN Interval
                Start		 8 Dec 2020 00:00:00.000000000
                Stop		 9 Dec 2020 00:00:00.000000000
            END Interval
            IntervalState		 Explicit
        END EVENTINTERVAL

        ConstConstraintsByStrands		 Yes
        UseSaveIntervalFile		 No
        SaveIntervalFile		 D:\Documents\Spacecraft-formation-control\numerical_studies\Letters_impuslive_control\Mission_design\STK\strand.int
        UseMinAngle		 No
        UseMaxAngle		 No
        UseMinLinkTime		 No
        LTDelayCriterion		 2
        TimeConvergence		 0.005
        AbsValueConvergence		 1e-14
        RelValueConvergence		 1e-08
        MaxTimeStep		 360
        MinTimeStep		 0.01
        UseLightTimeDelay		 Yes
        DetectEventsUsingSamplesOnly		 No
        UseLoadIntervalFile		 No
        BEGIN StrandObjIndexes
            StrandObj		 Place/Moscow
            StrandObj		 Satellite/Satellite2
        END StrandObjIndexes

        SaveMode		 1
        BEGIN StrandAccessesByIndex
            Strand		 0 1
            Start		  4.3927187648219733e+03
            Stop		  4.7339981239435629e+03
            Start		  1.0145455786918756e+04
            Stop		  1.0937299464369577e+04
            Start		  1.6014383194481288e+04
            Stop		  1.6856478948433152e+04
            Start		  2.1908512663225778e+04
            Stop		  2.2603356033293305e+04
            Start		  2.7811932700591529e+04
            Stop		  2.8208079596091837e+04
            Start		  3.3692756254600521e+04
            Stop		  3.3747811355302743e+04
            Start		  3.9230121508537231e+04
            Stop		  3.9631764951724748e+04
            Start		  4.4836743449370748e+04
            Stop		  4.5535112423735372e+04
            Start		  5.0586043696922898e+04
            Stop		  5.1429023641552332e+04
            Start		  5.6508244378765325e+04
            Stop		  5.7297304398265165e+04
            Start		  6.2720222995471675e+04
            Stop		  6.3044292252080442e+04
        END StrandAccessesByIndex


    END Definition

    BEGIN Extensions

        BEGIN ExternData
        END ExternData

        BEGIN ADFFileData
        END ADFFileData

        BEGIN Desc
            BEGIN ShortText

            END ShortText
            BEGIN LongText

            END LongText
        END Desc

        BEGIN Crdn
        END Crdn

        BEGIN Graphics

            BEGIN Attributes

                StaticColor		 #00ff00
                AnimationColor		 #00ffff
                AnimationLineWidth		 2
                StaticLineWidth		 3

            END Attributes

            BEGIN Graphics
                ShowGfx		 On
                ShowStatic		 Off
                ShowAnimationHighlight		 On
                ShowAnimationLine		 On
                ShowLinkDirection		 Off
            END Graphics
        END Graphics

        BEGIN VO
        END VO

    END Extensions

END Chain

